/***** This code was generated by Yaggo. Do not edit ******/

#ifndef __OMH_SKETCH_HPP__
#define __OMH_SKETCH_HPP__

#include <stdint.h>
#include <unistd.h>
#include <stdlib.h>
#include <getopt.h>
#include <errno.h>
#include <string.h>
#include <stdexcept>
#include <string>
#include <limits>
#include <vector>
#include <iostream>
#include <sstream>
#include <memory>

class omh_sketch {
 // Boiler plate stuff. Conversion from string to other formats
  static bool adjust_double_si_suffix(double &res, const char *suffix) { // adjust a double value based on a given SI (International System of Units) suffix
    if(*suffix == '\0') // empty string
      return true;
    if(*(suffix + 1) != '\0') // If the suffix contains more than one character, return False as the function only handle 1 character suffix
      return false;

    switch(*suffix) {
    case 'a': res *= 1e-18; break;
    case 'f': res *= 1e-15; break;
    case 'p': res *= 1e-12; break;
    case 'n': res *= 1e-9;  break;
    case 'u': res *= 1e-6;  break;
    case 'm': res *= 1e-3;  break;
    case 'k': res *= 1e3;   break;
    case 'M': res *= 1e6;   break;
    case 'G': res *= 1e9;   break;
    case 'T': res *= 1e12;  break;
    case 'P': res *= 1e15;  break;
    case 'E': res *= 1e18;  break;
    default: return false;
    }
    return true; // res is successfully adjusted
  }

  static double conv_double(const char *str, ::std::string &err, bool si_suffix) { // convert a string to a double value and handle possible errors
    char *endptr = 0; // initialize a character pointer endptr to null
    errno = 0; // errno is a global variable that stores the error number when a system call fails
    double res = strtod(str, &endptr); // converts the string pointed to by str to a double value. The conversion stops when the first invalid character is found. endptr is updated to point to this character
    if(endptr == str) { // no valid conversion could be performed
      err.assign("Invalid floating point string");
      return (double)0.0;
    }
    if(errno) { // not 0; means an error occurred during the conversion (for example, the number might be too large to represent)
      err.assign(strerror(errno)); // sets err to the error message corresponding to errno
      return (double)0.0;
    }
    bool invalid = // checks if there are any invalid characters in the string
      si_suffix ? !adjust_double_si_suffix(res, endptr) : *endptr != '\0'; // If this function returns false (indicating an invalid suffix), invalid is set to true
    if(invalid) {
      err.assign("Invalid character");
      return (double)0.0;
    }
    return res;
  }

  static int conv_enum(const char* str, ::std::string& err, const char* const strs[]) { // convert a string to an enumeration value; const char* const strs[]: This declares strs as an array of constant pointers to constant characters. Neither the pointers in the array nor the characters they point to can be changed
    int res = 0; // initializes an integer res to 0 (res keeps track of the current index in the array). This will be the result of the function if the string matches one of the strings in the array
    for(const char* const* cstr = strs; *cstr; ++cstr, ++res) // const char* const* cstr = strs;: This initializes a pointer cstr to the start of the array strs; *cstr;: This is the condition of the for loop. The loop will continue as long as the pointer that cstr points to is not null
      if(!strcmp(*cstr, str)) // If a match is found
        return res; // returns the current index res
    err += "Invalid constant '"; // constructs an error message if no match is found
    err += str;
    err += "'. Expected one of { ";
    for(const char* const* cstr = strs; *cstr; ++cstr) {
      if(cstr != strs)
        err += ", ";
      err += *cstr; // appending each string to err
    }
    err += " }";
    return -1; // The array of strings represents the possible values of the enumeration. If the string does not match any of the enumeration values, the function returns -1 and sets err to a descriptive error message
  }

  template<typename T>
  static bool adjust_int_si_suffix(T &res, const char *suffix) {
    if(*suffix == '\0') // empty string 
      return true;
    if(*(suffix + 1) != '\0') // contains more than one character
      return false;

    switch(*suffix) {
    case 'k': res *= (T)1000; break;
    case 'M': res *= (T)1000000; break;
    case 'G': res *= (T)1000000000; break;
    case 'T': res *= (T)1000000000000; break;
    case 'P': res *= (T)1000000000000000; break;
    case 'E': res *= (T)1000000000000000000; break;
    default: return false;
    }
    return true;
  }

  template<typename T>
  static T conv_int(const char *str, ::std::string &err, bool si_suffix) { // converts a string to an integer value of type T
    char *endptr = 0; // initialize a character pointer endptr to null
    errno = 0; // set the error number errno to 0
    long long int res = strtoll(str, &endptr, 0); // converts the string pointed to by str to a long long integer value; endptr is updated to point to this character
    if(endptr == str) { //  no valid conversion
      err.assign("Invalid signed int string");
      return (T)0;
    }
    if(errno) { // an error occurred during the conversion
      err.assign(strerror(errno)); // 
      return (T)0;
    }
    bool invalid =
      si_suffix ? !adjust_int_si_suffix(res, endptr) : *endptr != '\0'; // ì si_suffix, converts int to si suffix
    if(invalid) {
      err.assign("Invalid character");
      return (T)0;
    }
    if(res > ::std::numeric_limits<T>::max() || // checks if res is within the range of the type T
       res < ::std::numeric_limits<T>::min()) {
      err.assign("Value out of range");
      return (T)0;
    }
    return (T)res;
  }

  template<typename T>
  static T conv_uint(const char *str, ::std::string &err, bool si_suffix) { // converts a string to an unsigned integer value of type T
    char *endptr = 0; // initialize a character pointer endptr to null
    errno = 0; // set the error number errno to 0
    while(isspace(*str)) { ++str; } // skips any leading whitespace in the string
    if(*str == '-') { // checks if the string starts with a minus sign, which would indicate a negative number
      err.assign("Negative value");
      return (T)0;
    }
    unsigned long long int res = strtoull(str, &endptr, 0); // converts the string pointed to by str to an unsigned long long integer value. The conversion stops when the first invalid character is found. endptr is updated to point to this character
    if(endptr == str) { // no valid conversion
      err.assign("Invalid unsigned int string");
      return (T)0;
    }
    if(errno) { // error occurred during the conversion 
      err.assign(strerror(errno));
      return (T)0;
    }
    bool invalid =
      si_suffix ? !adjust_int_si_suffix(res, endptr) : *endptr != '\0';
    if(invalid) {
      err.assign("Invalid character");
      return (T)0;
    }
    if(res > ::std::numeric_limits<T>::max()) {
      err.assign("Value out of range");
      return (T)0;
    }
    return (T)res;
  }

  template<typename T>
  static ::std::string vec_str(const std::vector<T> &vec) { // converts a vector of elements of type T to a string
    ::std::ostringstream os; // creates an output string stream os
    for(typename ::std::vector<T>::const_iterator it = vec.begin(); // initializes a constant iterator it to the beginning of the vector vec
        it != vec.end(); ++it) {
      if(it != vec.begin()) // If the current element is not the first element in the vector, it inserts a comma into the stream os
        os << ",";
      os << *it; // inserts the current element into the stream os
    }
    return os.str(); // converts the contents of the stream os to a string and returns it
  }

  class string : public ::std::string { // derived class that inherits from the ::std::string class of the C++ Standard Library
  public:
    string() : ::std::string() {} // The default constructor string() which initializes the string to an empty string
    explicit string(const ::std::string &s) : std::string(s) {} // initializes the string with a copy of another ::std::string
    explicit string(const char *s) : ::std::string(s) {} // initializes the string with a C-string (null-terminated character array)
    int as_enum(const char* const strs[]) { // converts the string to an enumeration value
      ::std::string err;
      int res = conv_enum((const char*)this->c_str(), err, strs); // calling the c_str method on the current object (this), The c_str method is a member function of std::string that returns a pointer to a null-terminated character array with data equivalent to those stored in the string, converting the std::string to a C-string
      if(!err.empty())
        throw ::std::runtime_error(err);
      return res;
    }


    uint32_t as_uint32_suffix() const { return as_uint32(true); } // calls the as_uint32 method with si_suffix set to true, interpret any SI suffix in the string (like ‘k’ for kilo, ‘M’ for mega, etc.) when performing the conversion
    uint32_t as_uint32(bool si_suffix = false) const { // converts the string to a uint32_t
      ::std::string err;
      uint32_t res = conv_uint<uint32_t>((const char*)this->c_str(), err, si_suffix); // If conv_uint<uint32_t> returns an error message (i.e., the string cannot be converted to a uint32_t), the method constructs a detailed error message and throws a runtime error with this message
      if(!err.empty()) {
        ::std::string msg("Invalid conversion of '");
        msg += *this;
        msg += "' to uint32_t: ";
        msg += err;
        throw ::std::runtime_error(msg);
      }
      return res; // If no error occurred, the method returns the converted value
    }
    uint64_t as_uint64_suffix() const { return as_uint64(true); }
    uint64_t as_uint64(bool si_suffix = false) const {
      ::std::string err;
      uint64_t res = conv_uint<uint64_t>((const char*)this->c_str(), err, si_suffix);
      if(!err.empty()) {
        ::std::string msg("Invalid conversion of '");
        msg += *this;
        msg += "' to uint64_t: ";
        msg += err;
        throw ::std::runtime_error(msg);
      }
      return res;
    }
    int32_t as_int32_suffix() const { return as_int32(true); }
    int32_t as_int32(bool si_suffix = false) const {
      ::std::string err;
      int32_t res = conv_int<int32_t>((const char*)this->c_str(), err, si_suffix);
      if(!err.empty()) {
        ::std::string msg("Invalid conversion of '");
        msg += *this;
        msg += "' to int32_t: ";
        msg += err;
        throw ::std::runtime_error(msg);
      }
      return res;
    }
    int64_t as_int64_suffix() const { return as_int64(true); }
    int64_t as_int64(bool si_suffix = false) const {
      ::std::string err;
      int64_t res = conv_int<int64_t>((const char*)this->c_str(), err, si_suffix);
      if(!err.empty()) {
        ::std::string msg("Invalid conversion of '");
        msg += *this;
        msg += "' to int64_t: ";
        msg += err;
        throw ::std::runtime_error(msg);
      }
      return res;
    }
    int as_int_suffix() const { return as_int(true); }
    int as_int(bool si_suffix = false) const {
      ::std::string err;
      int res = conv_int<int>((const char*)this->c_str(), err, si_suffix);
      if(!err.empty()) {
        ::std::string msg("Invalid conversion of '");
        msg += *this;
        msg += "' to int_t: ";
        msg += err;
        throw ::std::runtime_error(msg);
      }
      return res;
    }
    long as_long_suffix() const { return as_long(true); }
    long as_long(bool si_suffix = false) const {
      ::std::string err;
      long res = conv_int<long>((const char*)this->c_str(), err, si_suffix);
      if(!err.empty()) {
        ::std::string msg("Invalid conversion of '");
        msg += *this;
        msg += "' to long_t: ";
        msg += err;
        throw ::std::runtime_error(msg);
      }
      return res;
    }
    double as_double_suffix() const { return as_double(true); }
    double as_double(bool si_suffix = false) const {
      ::std::string err;
      double res = conv_double((const char*)this->c_str(), err, si_suffix);
      if(!err.empty()) {
        ::std::string msg("Invalid conversion of '");
        msg += *this;
        msg += "' to double_t: ";
        msg += err;
        throw ::std::runtime_error(msg);
      }
      return res;
    }
  };

public:
  uint32_t                       k_arg;
  bool                           k_given;
  uint32_t                       l_arg;
  bool                           l_given;
  uint32_t                       m_arg;
  bool                           m_given;
  bool                           rc_flag;
  bool                           case_flag;
  const char *                   output_arg;
  bool                           output_given;
  const char *                   seedin_arg;
  bool                           seedin_given;
  const char *                   seedout_arg;
  bool                           seedout_given;
  bool                           randseed_flag;
  ::std::vector<const char *>    sequence_arg;
  typedef ::std::vector<const char *>::iterator sequence_arg_it;
  typedef ::std::vector<const char *>::const_iterator sequence_arg_const_it;

  enum { // defines six constants with values from 1000 to 1005
    START_OPT = 1000,
    RC_OPT,
    CASE_OPT,
    SEEDIN_OPT,
    SEEDOUT_OPT,
    RANDSEED_OPT
  };

  omh_sketch() :
    k_arg((uint32_t)5), k_given(false),
    l_arg((uint32_t)2), l_given(false),
    m_arg((uint32_t)1000), m_given(false),
    rc_flag(false),
    case_flag(false),
    output_arg(""), output_given(false),
    seedin_arg(""), seedin_given(false),
    seedout_arg(""), seedout_given(false),
    randseed_flag(false),
    sequence_arg() // default-initializes the sequence_arg member. Since sequence_arg is a vector, this means it is initialized to an empty vector
  { }

  omh_sketch(int argc, char* argv[]) :
    k_arg((uint32_t)5), k_given(false),
    l_arg((uint32_t)2), l_given(false),
    m_arg((uint32_t)1000), m_given(false),
    rc_flag(false),
    case_flag(false),
    output_arg(""), output_given(false),
    seedin_arg(""), seedin_given(false),
    seedout_arg(""), seedout_given(false),
    randseed_flag(false),
    sequence_arg()
  { parse(argc, argv); } // parse the command-line arguments and set the members of the class accordingly

  void parse(int argc, char* argv[]) {
    static struct option long_options[] = { // an array of option structures
      {"", 1, 0, 'k'}, // name: "", has_arg: 0=no_arg, 1=required_arg, 2=optional_arg; *flag=if not equal 0, getopt_long will set the value the fourth field to this variable; val = value to return, or to load into the variable pointed to by flag
      {"", 1, 0, 'l'},
      {"", 1, 0, 'm'},
      {"rc", 0, 0, RC_OPT},
      {"case", 0, 0, CASE_OPT},
      {"output", 1, 0, 'o'},
      {"seedin", 1, 0, SEEDIN_OPT},
      {"seedout", 1, 0, SEEDOUT_OPT},
      {"randseed", 0, 0, RANDSEED_OPT},
      {"help", 0, 0, 'h'},
      {"usage", 0, 0, 'U'},
      {"version", 0, 0, 'V'},
      {0, 0, 0, 0}
    };
    static const char *short_options = "hVUk:l:m:o:"; // short options (command-line options that start with -). Each character in the string is a short option. If a character is followed by a :, then the option requires an argument

    ::std::string err;
#define CHECK_ERR(type,val,which) if(!err.empty()) { ::std::cerr << "Invalid " #type " '" << val << "' for [" which "]: " << err << "\n"; exit(1); } // prints an error message to the standard error stream (std::cerr). The #type is a stringizing operator that converts the macro parameter type to a string literal. The error message includes the type, the val, the which, and the error string err
    while(true) { // uses the getopt_long function to parse command-line arguments
      int index = -1;
      // returns the option character for a short option or the value in the val field for a long option. If getopt_long encounters a long option, it sets index to the index of the option in the long_options array
      int c = getopt_long(argc, argv, short_options, long_options, &index); // initializes an integer index to -1. This variable will be used to store the index of the long option in the long_options array when getopt_long encounters a long option
      if(c == -1) break; // checks if getopt_long has processed all command-line arguments. If it has, getopt_long returns -1 and the loop is exited
      switch(c) {
      case ':':
        ::std::cerr << "Missing required argument for " // prints an error message to the standard error stream and exits the program with an error status code. The error message includes the name of the option that is missing an argument
                  << (index == -1 ? ::std::string(1, (char)optopt) : std::string(long_options[index].name))
                  << ::std::endl;
        exit(1);
      case 'h':
        ::std::cout << usage() << "\n\n" << help() << std::endl; // prints the usage message and the help message, and exits the program with status code 0
        exit(0);
      case 'U':
        ::std::cout << usage() << "\nUse --help for more information." << std::endl; // prints the usage message and exits the program with status code 0
        exit(0);
      case 'V':
        print_version(); // prints the version information and exits the program with status code 0
        exit(0);
      case '?':
        ::std::cerr << "Use --usage or --help for some help\n"; // unknown option is encountered, prints a help message and exits the program with status code 1
        exit(1);
      case 'k':
        k_given = true;
        k_arg = conv_uint<uint32_t>((const char*)optarg, err, false); // convert the option argument to an unsigned 32-bit integer and store it in the corresponding member variable
        CHECK_ERR(uint32_t, optarg, "-k")
        break;
      case 'l':
        l_given = true;
        l_arg = conv_uint<uint32_t>((const char*)optarg, err, false);
        CHECK_ERR(uint32_t, optarg, "-l")
        break;
      case 'm':
        m_given = true;
        m_arg = conv_uint<uint32_t>((const char*)optarg, err, false);
        CHECK_ERR(uint32_t, optarg, "-m")
        break;
      case RC_OPT:
        rc_flag = true;
        break;
      case CASE_OPT:
        case_flag = true;
        break;
      case 'o':
        output_given = true;
        output_arg = optarg;
        break;
      case SEEDIN_OPT:
        seedin_given = true;
        seedin_arg = optarg;
        break;
      case SEEDOUT_OPT:
        seedout_given = true;
        seedout_arg = optarg;
        break;
      case RANDSEED_OPT:
        randseed_flag = true;
        break;
      }
    }

    // Parse arguments
    if(argc - optind < 0) // checks if the number of command-line arguments (argc) minus the current option index (optind) is less than 0; The optind variable is from the getopt library and it’s the index of the next element of the argv[] array to be processed; If this condition is true, it means there are fewer command-line arguments than expected
      error("Requires at least 0 argument."); // calls the error function with a string argument to print an error message and likely exit the program
    for( ; optind < argc; ++optind) {
      sequence_arg.push_back(argv[optind]);
    }
  }
  static const char * usage() { return "Usage: omh_sketch [options] sequence:PATH+"; }
  class error {
    int code_;
    std::ostringstream msg_;

    // Select the correct version (GNU or XSI) version of
    // strerror_r. strerror_ behaves like the GNU version of strerror_r,
    // regardless of which version is provided by the system.
    static const char* strerror__(char* buf, int res) {
      return res != -1 ? buf : "Invalid error";
    }
    static const char* strerror__(char* buf, char* res) {
      return res;
    }
    static const char* strerror_(int err, char* buf, size_t buflen) {
      return strerror__(buf, strerror_r(err, buf, buflen));
    }
    struct no_t { };

  public:
    static no_t no;
    error(int code = EXIT_FAILURE) : code_(code) { }
    explicit error(const char* msg, int code = EXIT_FAILURE) : code_(code)
      { msg_ << msg; }
    error(const std::string& msg, int code = EXIT_FAILURE) : code_(code)
      { msg_ << msg; }
    error& operator<<(no_t) {
      char buf[1024];
      msg_ << ": " << strerror_(errno, buf, sizeof(buf));
      return *this;
    }
    template<typename T>
    error& operator<<(const T& x) { msg_ << x; return (*this); }
    ~error() {
      ::std::cerr << "Error: " << msg_.str() << "\n"
                  << usage() << "\n"
                  << "Use --help for more information"
                  << ::std::endl;
      exit(code_);
    }
  };
  static const char * help() { return
    "Compute sketches of sequence in fasta file\n\n"
    "Options (default value in (), *required):\n"
    " -k                                       kmer length (5)\n"
    " -l                                       vector length (2)\n"
    " -m                                       sketch length (1000)\n"
    "     --rc                                 include reverse complement (false)\n"
    "     --case                               do not ignore case (false)\n"
    " -o, --output=PATH                        Write sketches to given file instead of 1 file per fasta entry\n"
    "     --seedin=PATH                        Seed in\n"
    "     --seedout=PATH                       Seed out\n"
    "     --randseed                           If no seedin, generate a random seed (defaults to default seed) (false)\n"
    " -U, --usage                              Usage\n"
    " -h, --help                               This message\n"
    " -V, --version                            Version";
  }
  static const char* hidden() { return ""; }
  void print_version(::std::ostream &os = std::cout) const {
#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.0.0"
#endif
    os << PACKAGE_VERSION << "\n";
  }
  void dump(::std::ostream &os = std::cout) {
    os << "k_given:" << k_given << " k_arg:" << k_arg << "\n";
    os << "l_given:" << l_given << " l_arg:" << l_arg << "\n";
    os << "m_given:" << m_given << " m_arg:" << m_arg << "\n";
    os << "rc_flag:" << rc_flag << "\n";
    os << "case_flag:" << case_flag << "\n";
    os << "output_given:" << output_given << " output_arg:" << output_arg << "\n";
    os << "seedin_given:" << seedin_given << " seedin_arg:" << seedin_arg << "\n";
    os << "seedout_given:" << seedout_given << " seedout_arg:" << seedout_arg << "\n";
    os << "randseed_flag:" << randseed_flag << "\n";
    os << "sequence_arg:" << vec_str(sequence_arg) << "\n";
  }
};
#endif // __OMH_SKETCH_HPP__"
