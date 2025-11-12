#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <array>
#include <cstdint>
#include <random>
#include <string>
#include <vector>
#include <algorithm>
#include <stdexcept>

#include <boost/dynamic_bitset.hpp>

#include <mpi.h>

// SQD-HPC
#include <qiskit/addon/sqd/subsampling.hpp>
#include <qiskit/addon/sqd/configuration_recovery.hpp>
#include <qiskit/addon/sqd/fermion.hpp>

// SBD
#include "sbd/sbd.h"

namespace py = pybind11;

static std::vector<boost::dynamic_bitset<>> to_bitsets(
    const std::vector<std::string> &bitstrings
) {
    std::vector<boost::dynamic_bitset<>> out;
    out.reserve(bitstrings.size());
    for (const auto &s : bitstrings) {
        out.emplace_back(s);
    }
    return out;
}

static std::string bitset_to_string(const boost::dynamic_bitset<> &bs) {
    std::string s;
    s.reserve(bs.size());
    for (std::size_t i = 0; i < bs.size(); ++i) {
        s.push_back(bs.test(i) ? '1' : '0');
    }
    std::reverse(s.begin(), s.end());
    return s;
}

static std::pair<std::vector<double>, std::vector<double>>
compute_avg_occs_from_samples(
    const std::vector<boost::dynamic_bitset<>> &bitstrings,
    const std::vector<double> &weights
) {
    if (bitstrings.empty()) return {{}, {}};
    const std::size_t nbits = bitstrings[0].size();
    if (nbits % 2 != 0) {
        throw std::runtime_error("bitstrings must have even length (beta||alpha)");
    }
    const std::size_t norb = nbits / 2;
    std::vector<double> alpha(norb, 0.0), beta(norb, 0.0);
    for (std::size_t k = 0; k < bitstrings.size(); ++k) {
        const auto &bs = bitstrings[k];
        double p = weights[k];
        for (std::size_t i = 0; i < norb; ++i) {
            if (bs.test(i))        beta[i]  += p;
            if (bs.test(i + norb)) alpha[i] += p;
        }
    }
    return {alpha, beta};
}

py::list run_sqd_only(
    const std::vector<std::string> &bitstrings,
    const std::vector<double> &weights,
    const std::vector<double> &avg_occs_alpha,
    const std::vector<double> &avg_occs_beta,
    std::uint64_t n_alpha,
    std::uint64_t n_beta,
    std::size_t samples_per_batch = 1000,
    std::size_t num_batches       = 1,
    std::uint64_t max_ci_dim      = 0,
    std::uint64_t seed            = 12345
) {
    if (bitstrings.size() != weights.size()) {
        throw std::runtime_error("bitstrings and weights must have the same length");
    }

    auto bits = to_bitsets(bitstrings);

    std::vector<double> aocc = avg_occs_alpha;
    std::vector<double> bocc = avg_occs_beta;
    if (aocc.empty() || bocc.empty()) {
        auto tmp = compute_avg_occs_from_samples(bits, weights);
        aocc = std::move(tmp.first);
        bocc = std::move(tmp.second);
    }

    std::array<std::vector<double>, 2> avg_occs{aocc, bocc};
    std::array<std::uint64_t, 2> nelec{n_alpha, n_beta};

    std::mt19937_64 rng(seed);

    auto [rec_bits, rec_wts] = Qiskit::addon::sqd::recover_configurations(
        bits,
        weights,
        avg_occs,
        nelec,
        rng
    );

    std::size_t available = rec_bits.size();
    unsigned int spb = static_cast<unsigned int>(samples_per_batch);
    if (available < spb) spb = static_cast<unsigned int>(available);
    if (spb == 0) {
        throw std::runtime_error("SQD: no configurations available for subsampling.");
    }

    auto batches = Qiskit::addon::sqd::subsample_multiple_batches(
        rec_bits,
        rec_wts,
        spb,
        static_cast<unsigned int>(num_batches),
        rng
    );

    py::list py_batches;
    for (const auto &batch : batches) {
        std::optional<unsigned int> max_dim_opt = std::nullopt;
        if (max_ci_dim > 0) max_dim_opt = static_cast<unsigned int>(max_ci_dim);

        auto ci_bitsets =
            Qiskit::addon::sqd::bitstrings_to_ci_strings_symmetrize_spin(
                batch,
                max_dim_opt,
                std::nullopt
            );

        py::list py_ci;
        for (const auto &ci_bs : ci_bitsets) {
            py_ci.append(bitset_to_string(ci_bs));
        }
        py_batches.append(py_ci);
    }

    return py_batches;
}

py::dict run_sqd_sbd(
    const std::vector<std::string> &bitstrings,
    const std::vector<double> &weights,
    std::uint64_t n_alpha,
    std::uint64_t n_beta,
    const std::string &fcidump_path,
    const std::vector<double> &avg_occs_alpha = {},
    const std::vector<double> &avg_occs_beta  = {},
    std::size_t samples_per_batch             = 1000,
    std::size_t num_batches                   = 1,
    std::uint64_t /*max_ci_dim*/              = 0,
    std::size_t bit_length                    = sizeof(std::size_t)*8,
    std::uint64_t seed                        = 12345,
    bool beta_first                           = true   
) {
    int mpi_initialized = 0;
    MPI_Initialized(&mpi_initialized);
    if (!mpi_initialized) {
        MPI_Init(nullptr, nullptr);
    }

    auto bs_measured = to_bitsets(bitstrings); // strings to bitsets

    std::vector<double> aocc = avg_occs_alpha;
    std::vector<double> bocc = avg_occs_beta;
    if (aocc.empty() || bocc.empty()) {
        auto est = compute_avg_occs_from_samples(bs_measured, weights);
        aocc = std::move(est.first);
        bocc = std::move(est.second);
    }

    std::array<std::vector<double>, 2> avg_occs{aocc, bocc};
    std::array<std::uint64_t, 2> nelec{n_alpha, n_beta};

    std::mt19937_64 rng(seed);

    // sqd recovery
    auto [rec_bits, rec_wts] = Qiskit::addon::sqd::recover_configurations(
        bs_measured,
        weights,
        avg_occs,
        nelec,
        rng
    );
    if (rec_bits.empty()) {
        throw std::runtime_error("SQD: no configurations after recovery.");
    }

    const std::size_t full_len = rec_bits[0].size();
    if (full_len % 2 != 0) {
        throw std::runtime_error("Recovered bitstrings must have even length.");
    }
    const std::size_t L = full_len / 2;

    // collect unique alpha and beta determinants
    std::vector<std::string> alpha_strings;
    std::vector<std::string> beta_strings;

    for (const auto &bs : rec_bits) {
        std::string full;
        full.reserve(full_len);
        for (std::size_t i = 0; i < full_len; ++i) {
            full.push_back(bs.test(i) ? '1' : '0');
        }
        std::reverse(full.begin(), full.end());

        std::string left  = full.substr(0, L);
        std::string right = full.substr(L);

        std::string beta  = beta_first ? left  : right;
        std::string alpha = beta_first ? right : left;

        alpha_strings.push_back(std::move(alpha));
        beta_strings.push_back(std::move(beta));
    }

    std::sort(alpha_strings.begin(), alpha_strings.end());
    alpha_strings.erase(std::unique(alpha_strings.begin(), alpha_strings.end()), alpha_strings.end());

    std::sort(beta_strings.begin(), beta_strings.end());
    beta_strings.erase(std::unique(beta_strings.begin(), beta_strings.end()), beta_strings.end());

    // turn into sbd format
    std::vector<std::vector<size_t>> adet;
    std::vector<std::vector<size_t>> bdet;
    adet.reserve(alpha_strings.size());
    bdet.reserve(beta_strings.size());

    for (const auto &a : alpha_strings) {
        adet.push_back(sbd::from_string(a, bit_length, L));
    }
    for (const auto &b : beta_strings) {
        bdet.push_back(sbd::from_string(b, bit_length, L));
    }

    sbd::FCIDump fcidump = sbd::LoadFCIDump(fcidump_path);

    // diagonalizing with sbd
    sbd::tpb::SBD sbd_data;
    double energy = 0.0;
    std::vector<double> density;
    std::vector<std::vector<size_t>> carryover_bitstrings;
    std::vector<std::vector<double>> one_p_rdm;
    std::vector<std::vector<double>> two_p_rdm;

    sbd::tpb::diag(
        MPI_COMM_WORLD,
        sbd_data,
        fcidump,
        adet,
        bdet,
        "",
        "",
        energy,
        density,
        carryover_bitstrings,
        one_p_rdm,
        two_p_rdm
    );

    py::dict out;
    out["energy"] = energy;
    out["density"] = density;

    py::list carryover_py;
    for (const auto &cbs : carryover_bitstrings) {
        std::string s = sbd::makestring(cbs, bit_length, L);
        carryover_py.append(s);
    }
    out["carryover_bitstrings"] = carryover_py;
    return out;
}

void init_sqd_sbd(py::module_ &m) {
    m.def(
        "run_sqd_only",
        &run_sqd_only,
        py::arg("bitstrings"),
        py::arg("weights"),
        py::arg("avg_occs_alpha") = std::vector<double>{},
        py::arg("avg_occs_beta")  = std::vector<double>{},
        py::arg("n_alpha"),
        py::arg("n_beta"),
        py::arg("samples_per_batch") = 1000,
        py::arg("num_batches")       = 1,
        py::arg("max_ci_dim")        = 0,
        py::arg("seed")              = 12345
    );

    m.def(
        "run_sqd_sbd",
        &run_sqd_sbd,
        py::arg("bitstrings"),
        py::arg("weights"),
        py::arg("n_alpha"),
        py::arg("n_beta"),
        py::arg("fcidump_path"),
        py::arg("avg_occs_alpha")    = std::vector<double>{},
        py::arg("avg_occs_beta")     = std::vector<double>{},
        py::arg("samples_per_batch") = 1000,
        py::arg("num_batches")       = 1,
        py::arg("max_ci_dim")        = 0,
        py::arg("bit_length")        = sizeof(std::size_t)*8,
        py::arg("seed")              = 12345,
        py::arg("beta_first")        = true 
    );
}

