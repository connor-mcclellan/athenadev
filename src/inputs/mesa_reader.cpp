#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>  // runtime_error

#include "../athena.hpp"

void MesaReader(const char *filename, std::vector<std::string> &headers,
                std::vector<std::string> &hdata, std::vector<std::string> &vars,
                std::vector<double> &data, bool verbose = false) {
  // \brief Read a single 1D MESA profile data dump

  std::ifstream file(filename, std::ios::in);
  std::string line;
  std::stringstream stream;

  // Parse the header data as strings (header is mixed type)
  // -------------------------------------------------------

  // Count number of header entries
  // This is done to avoid reallocation as new elements are added - instead,
  // count once, allocate once, and then fill the vector
  std::getline(file, line);
  stream.str(line);
  std::string substr;
  int nhead = 0;
  while (stream) {
    substr.clear();
    stream >> substr;
    if (substr[0] != '\0') {
      nhead++;
      if (verbose)
        std::cout << substr << std::endl;
    }
  }

  // Resize string vectors to contain the header info
  headers.reserve(nhead);
  hdata.reserve(nhead);

  // Parse header names
  std::getline(file, line);
  stream.clear();
  stream.str(line);
  int h = 0;
  while (stream) {
    substr.clear();
    stream >> substr;
    if (substr[0] != '\0') {
      headers.push_back(substr);
      h++;
      if (verbose)
        std::cout << substr << std::endl;
    }
  }

  // Parse header values
  std::getline(file, line);
  stream.clear();
  stream.str(line);
  h = 0;
  while (stream) {
    substr.clear();
    stream >> substr;
    if (substr[0] != '\0') {
      hdata.push_back(substr);
      h++;
      if (verbose)
        std::cout << substr << std::endl;
    }
  }

  // Skip blank lines after header
  while ((std::getline(file, line)) && (line.empty()))
    continue;

  // Count number of variable data columns
  stream.clear();
  stream.str(line);
  int nvar = 0;
  while (stream) {
    substr.clear();
    stream >> substr;
    if (substr[0] != '\0') {
      nvar++;
      if (verbose)
        std::cout << substr << std::endl;
    }
  }

  // Resize vector for variable names
  vars.reserve(nvar);

  // Parse variable names
  std::getline(file, line);
  stream.clear();
  stream.str(line);
  int n = 0;
  while (stream) {
    substr.clear();
    stream >> substr;
    if (substr[0] != '\0') {
      if (verbose)
        std::cout << substr << std::endl;
      vars.push_back(substr);
      n++;
    }
  }

  // Parse the variable data - m positions x n variables
  // ---------------------------------------------------

  // Count number of rows in the data block
  int npos = 1;
  while (std::getline(file, line) && (file.peek() != EOF)) {
    npos++;
  }

  // Return to beginning of file
  file.clear();            // clear the eof flag so we can still read the file
  file.seekg(0, file.beg); // reset position to beginning of file

  // Skip to the blank line separating the header and the data
  while ((std::getline(file, line)) && (!line.empty()))
    continue;

  // Skip the rows containing the variable index and variable name
  std::getline(file, line);
  std::getline(file, line);

  // Now parse the data block
  data.resize(npos * nvar, 0.0); // resize is needed rather than reserve, so that the vector contains a nonzero number of items
  stream.precision(17);
  for (int m = 0; m < npos; m++) {
    std::getline(file, line);
    stream.clear();
    stream.str(line);
    for (int n = 0; n < nvar; n++) {
      substr.clear();
      stream >> substr;
      data.at(m + npos * (n)) =
          std::stod(substr); // Note the index ordering - fastest
                             // access should be to position data
      if (verbose)
        std::cout << "n: " << n << "    m: " << m << "    var: " << substr
                  << std::endl;
    }
  }
  file.close();
}

void GetMesaData(std::string field, std::vector<std::string> vars,
                 std::vector<double> vdata, double *pdata) {
  // \brief Copies M=(size of `vdata`)/(size of `vars`) elements from `vdata`
  // into `*pdata`, starting at the index n at which `field` and `vars[n]` match
  int nradii = vdata.size() / vars.size();
  for (int n = 0; n < vars.size(); n++) {
    if (field.compare(vars[n]) == 0) {
      for (int m = 0; m < nradii; m++) {
        pdata[m] = vdata[m + nradii * (n)];
      }
      return;
    }
  }
  std::stringstream msg;
  msg << "### FATAL ERROR in GetMesaData" << std::endl
      << "Field '" << field << "' not found in MESA profile." << std::endl;
  ATHENA_ERROR(msg);
  return;
}

void GetMesaHeader(std::string field, std::vector<std::string> headers,
                   std::vector<std::string> hdata, double &value) {
  // \brief Converts a string from a MESA header field into a double
  for (int h = 0; h < headers.size(); h++) {
    if (field.compare(headers[h]) == 0) {
      value = std::stod(hdata[h]);
      return;
    }
  }
  std::stringstream msg;
  msg << "### FATAL ERROR in GetMesaData" << std::endl
      << "Field '" << field << "' not found in MESA profile." << std::endl;
  ATHENA_ERROR(msg);
  return;
}
