#ifndef MESA_READER_HPP_
#define MESA_READER_HPP_

#include <string>
#include <vector>
#include <stdexcept>  // runtime_error

#include "../athena.hpp"             // Real
#include "../athena_arrays.hpp"      // AthenaArray
#include "../utils/interp_table.hpp" // InterpTable2D

void MesaReader(const char *filename, std::vector<std::string> &headers,
                std::vector<std::string> &hdata, std::vector<std::string> &vars,
                std::vector<double> &data, bool verbose = false);

// Functions for retrieving MESA data fields by name
void GetMesaData(std::string field, std::vector<std::string> vars,
                 std::vector<Real> vdata, Real *pdata);
void GetMesaHeader(std::string field, std::vector<std::string> headers,
                   std::vector<std::string> hdata, Real &value);


#endif
