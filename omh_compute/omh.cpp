
#include <omh.hpp>
#include <iostream>

char comp(char c) {
  switch(c) {
  case 'a': return 't';
  case 'c': return 'g';
  case 'g': return 'c';
  case 't': return 'a';
  case 'A': return 'T';
  case 'C': return 'G';
  case 'G': return 'C';
  case 'T': return 'A';
  }
  return c;
}

std::string reverse_complement(const std::string& sequence) { // return reverse complement of a sequence
  std::string res(sequence.size(), '\0'); // creates a string res of the same size as sequence, initially filled with null characters

  auto sf = sequence.begin(); // create iterators pointing to the start and end of sequence and res
  auto sl = sequence.end() - 1;
  auto rf = res.begin();
  auto rl = res.end() - 1;

  for( ; sf < sl; ++sf, --sl, ++rf, --rl) { // iterates over sequence and res from both ends towards the middle. For each pair of elements, it computes the complement of the element from sequence and stores it in the corresponding position in res
    *rf = comp(*sl);
    *rl = comp(*sf);
  }
  if(sf == sl) {
    *rf = comp(*sf);
  }
  return res;
}

double compare_sketch_pair(const char* p1, const char* p2, unsigned m, unsigned k, unsigned l, bool circular) { // compares two sketches of DNA sequences and returns a similarity score
  const unsigned block = std::max(l, (unsigned)1) * k; // computes the size of each sketch in the character arrays
  unsigned count = 0; // initializes a counter that will keep track of the number of matching sketches
  if(!circular || l < 2) { 
    for(unsigned i = 0; i < m; ++i, p1 += block, p2 += block)
      count += memcmp(p1, p2, block) == 0; // compare the sketches. If two sketches are identical, memcmp returns 0 and the counter is incremented
  } else {
    for(unsigned i = 0; i < m; ++i, p1 += block, p2 += block) { // compares the sketches in a circular manner, meaning that it considers rotations of the sketches as well. If a rotation of a sketch from p1 matches a sketch from p2, the counter is incremented
      for(unsigned j = 0; j < l; ++j) {
        if(memcmp(p1, p2 + j * k, block - j * k) == 0 && memcmp(p1 + block - j * k, p2, j * k) == 0) {
          ++count;
          break;
        }
      }
    }
  }
  return (double)count / m; // computes the similarity score as the number of matching sketches divided by the total number of sketches
}

double compare_sketches(const sketch& sk1, const sketch& sk2, ssize_t m, bool circular) { // ssize_t m: This is the number of sketches to compare
  if(sk1.k != sk2.k || sk1.l != sk2.l) return -1; // Different k or l
  if(m < 0) m = std::min(sk1.m, sk2.m); 
  if(m > sk1.m || m > sk2.m) return -1;  // Too short

  const unsigned block = std::max(sk1.l, (unsigned)1) * sk1.k;
  if(sk1.data.size() < m * block || sk2.data.size() < m * block) return -1; // Truncated
  const double fwd_score = compare_sketch_pair(sk1.data.data(), sk2.data.data(), m, sk1.k, sk1.l, circular); // computes a similarity score for the forward strands of the sequences using the compare_sketch_pair function

  double bwd_score = 0.0;
  if(!sk1.rcdata.empty()) {
    bwd_score = compare_sketch_pair(sk1.rcdata.data(), sk2.data.data(), m, sk1.k, sk1.l, circular); // computes a similarity score for the reverse complements of the sequences, if they are available
  } else if(!sk2.rcdata.empty()) {
    bwd_score = compare_sketch_pair(sk1.data.data(), sk2.rcdata.data(), m, sk1.k, sk1.l, circular); // If only one sequence has a reverse complement, it compares the forward strand of the other sequence with the reverse complement
  }

  return std::max(fwd_score, bwd_score); // the maximum of the forward and reverse complement similarity scores
}

void sketch::read(std::istream& is) { //  reads a sketch from an input stream and stores it in the sketch object; the sketch object contains the sketch read from the input stream
  if(is.get() != '>') throw std::runtime_error("Invalid format, expecting '>'"); // reads the first character from the input stream and checks if it’s a > character, which is a common convention in bioinformatics file formats. If the first character is not >, it throws an exception
  bool rc;
  is >> name >> k >> l >> m >> rc; // reads the name of the sequence and the parameters k, l, m, and rc from the input stream
  is.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // ignores the rest of the current line in the input stream, up to and including the newline character
  if(!is.good()) throw std::runtime_error("Invalid header"); 
  data.resize(std::max(l, (unsigned)1) * m * k); //  resizes the data vector to hold the sketch data
  is.read(data.data(), data.size()); // reads the sketch data from the input stream into the data vector
  is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  if(!is.good()) throw std::runtime_error("Error reading sketch");
  if(rc) {
    rcdata.resize(std::max(l, (unsigned)1) * m * k); 
    is.read(rcdata.data(), rcdata.size()); // reads a sketch of the reverse complement of the sequence from the input stream into the rcdata vector
    is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    if(!is.good()) throw std::runtime_error("Error reading rc sketch");
  }
}

void sketch::write(std::ostream& os) const { // writes a sketch to an output stream
  os << '>' << name << ' ' << k << ' ' << l << ' ' << m << ' ' << !rcdata.empty() << '\n'; // writes the name of the sequence and the parameters k, l, m, and a boolean indicating whether rcdata is not empty (which would mean that the sketch includes the reverse complement of the sequence) to the output stream
  os.write(data.data(), data.size()); // writes the sketch data to the output stream
  os << '\n';
  if(!rcdata.empty()) {
    os.write(rcdata.data(), rcdata.size()); // writes the sketch data for the reverse complement of the sequence to the output stream
    os << '\n';
  }
}
