#ifndef __OMH_H__
#define __OMH_H__

#include <cstring>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <random>
#include <limits>

#include <xxhash.hpp>
#include <seeded_prg.hpp>

std::string reverse_complement(const std::string&);

struct mer_info {
  """ mer_info class, initialize by mer_info(p, o, h)
  """
  size_t pos; // position in the sequence
  uint64_t hash; // store the hashing value
  unsigned occ; // store the number of occurences
  mer_info(size_t p, unsigned o, uint64_t h) // constructor
    : pos(p)
    , hash(h)
    , occ(o)
  { }
};

// Compute the position in sequence of the k-mers picked by omh, and
// passed them 1 by 1 to block. block takes 3 arguments: i \in [m], j
// \in [l] and the position.
template<typename EngineT, typename BT>
void omh_pos(const std::string& seq, unsigned k, unsigned l, unsigned m, EngineT& prg, BT block) {
  if(seq.size() < k) return;
  const bool weight = l > 0;
  if(l == 0) l = 1;

  std::vector<mer_info> mers;
  std::unordered_map<std::string, unsigned> occurrences; // default dictionary: counting occurence of k-mer
  size_t pos[l];

  //  Create list of k-mers with occurrence numbers
  for(size_t i = 0; i < seq.size() - k + 1; ++i) {
    auto occ = occurrences[seq.substr(i, k)]++;
    mers.emplace_back(i, occ, (uint64_t)0); // (pos_in_seq, occ, hash)
  }

  xxhash hash; // creates an instance of the xxhash class, which is used to compute hash values
  for(unsigned i = 0; i < m; ++i) {
    const auto seed = prg(); // prg: pseudo random generator
    for(auto& meri : mers) {
      hash.reset(seed); // resets the xxhash object with the new seed
      hash.update(&seq.data()[meri.pos], k); // updates the hash value with the k-mer starting at position meri.pos in the sequence. The seq.data() function returns a pointer to the sequence string, and &seq.data()[meri.pos] gives a pointer to the start of the k-mer.
      if(weight) hash.update(&meri.occ, sizeof(meri.occ)); // updates the hash value with the number of occurrences of the k-mer. The sizeof(meri.occ) function returns the size of the occ member in bytes.
      meri.hash = hash.digest(); // computes the final hash value and stores it in the hash member of the mer_info object
    }

    std::partial_sort(mers.begin(), mers.begin() + l, mers.end(), [&](const mer_info& x, const mer_info& y) { return x.hash < y.hash; }); // partially sorts the mers vector. The elements from mers.begin() to mers.begin() + l are sorted in ascending order based on their hash values. The rest of the vector is not sorted.
    std::sort(mers.begin(), mers.begin() + l, [&](const mer_info& x, const mer_info& y) { return x.pos < y.pos; }); // sorts the first l elements of the mers vector in ascending order based on their positions in the sequence.  the first l elements of mers are the l k-mers with the smallest hash values, and they are ordered by their positions in the sequence.
    for(unsigned j = 0; j < l; ++j)
      block(i, j, mers[j].pos);
  }
}

struct sketch { // a data structure that represents a sketch of a DNA sequence
  std::string       name; // holds the name of the sequence
  unsigned          k, l, m; // k: kmer-length, l: vector length, m: sketch length
  std::vector<char> data; // a vector of characters that likely holds the sketch data for the original sequence.
  std::vector<char> rcdata; //  a vector of characters that likely holds the sketch data for the reverse complement of the sequence

  bool operator==(const sketch& rhs) const { // an equality operator that checks if two sketches are equal. Two sketches are considered equal if their k, l, m, data, and rcdata members are all equal
    return k == rhs.k && l == rhs.l && m == rhs.m && data == rhs.data && rcdata == rhs.rcdata;
  }

  // Read a sketch into this
  void read(std::istream& is); // reads a sketch from an input stream and stores it in the sketch object

  // Read a sketch
  static sketch from_stream(std::istream& is) { // creates a new sketch object, reads a sketch from an input stream into it, and returns it
    sketch sk;
    sk.read(is);
    return sk;
  }

  // Write sketch
  void write(std::ostream& os) const; //  writes the sketch to an output stream.
};

template<typename EngineT = std::mt19937_64>
class omh_sketcher {
protected:
  unsigned           m_k, m_l, m_m;
  LongSeed<EngineT>& m_seed;
  EngineT            m_prg;

  inline void write_sketch(std::ostream& os, const std::string& seq) {
    omh_pos(seq, m_k, m_l, m_m, m_prg, [&os, &seq, this](unsigned i, unsigned j, size_t pos) { os.write(seq.data() + pos, m_k); }); // writes the k-mer starting at position pos in the sequence to the output stream
  }

  inline void compute_sketch(char* ptr, const char* seq) {
    omh_pos(seq, m_k, m_l, m_m, m_prg,
            [&ptr, &seq, this](unsigned i, unsigned j, size_t pos) { memcpy(ptr, seq + pos, m_k); ptr += m_k; }); //  the character array pointed to by ptr contains the sketches of all k-mers in the sequence.
  }

public:
  omh_sketcher(unsigned k, unsigned l, unsigned m, LongSeed<EngineT>& seed)
    : m_k(k), m_l(l), m_m(m)
    , m_seed(seed)
  { }

  // Compute the positions in seq of the k-mers picked by omh and append them to res
  void pos(const std::string& seq, std::vector<size_t>& res) {
    m_seed(m_prg);
    omh_pos(seq, m_k, m_l, m_m, m_prg, [&res](unsigned i, unsigned j, size_t pos) { res.push_back(pos); }); // res vector contains the positions of all k-mers in the sequence that were selected by the omh_pos function
  }

  // Compute the positions in seq and return vector
  inline std::vector<size_t> pos(const std::string& seq) {
    std::vector<size_t> res;
    pos(seq, res); // call pos above
    return res;
  }

  // Write sketch to an ostream
  void write(std::ostream& os, const std::string& name, const std::string& seq, bool rc = false) { // os contains a sketch of the sequence and possibly its reverse complement, along with some metadata
    m_seed(m_prg);
    os << '>' << name << ' ' << m_k << ' ' << m_l << ' ' << m_m << ' ' << rc << '\n';
    write_sketch(os, seq);
    os << '\n';
    if(rc) {
      m_seed(m_prg);
      std::string rcseq = reverse_complement(seq);
      write_sketch(os, rcseq);
      os << '\n';
    }
  }

  // Sketch from sequence
  void compute(const std::string& seq, sketch& sk, bool rc = false) {
    sk.k = m_k;
    sk.l = m_l;
    sk.m = m_m;
    sk.data.resize(std::max(sk.l, (unsigned)1) * sk.m * sk.k);
    m_seed(m_prg);
    compute_sketch(sk.data.data(), seq.data()); // computes a sketch of a DNA sequence and stores it in a sketch object. The sketch object is updated with the parameters k, l, m, and the sketch data.

    if(rc) {
      m_seed(m_prg);
      std::string rcseq = reverse_complement(seq);
      sk.rcdata.resize(sk.l * sk.m * sk.k);
      compute_sketch(sk.rcdata.data(), rcseq.data());
    }
  }

  // Sketch from sequence
  sketch compute(const std::string& seq, bool rc = false) {
    sketch sk;
    compute(seq, sk, rc);
    return sk;
  }

  double compare(const std::string& s1, const std::string& s2, bool rc = false) {
    const auto sk1 = compute(s1, rc);
    const auto sk2 = compute(s2, rc); // computes sketches of two DNA sequences (and their reverse complements if rc is true)
    return compare_sketches(sk1, sk2); // returns a double value that represents the similarity between the two sequences
  }
};


// Compare 2 sketches. Return agreement in [0, 1]
double compare_sketches(const sketch& sk1, const sketch& sk2, ssize_t m = -1, bool circular = false);


#endif /* __OMH_H__ */
