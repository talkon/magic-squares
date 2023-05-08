#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

typedef std::vector<size_t> vec_t;
typedef std::vector<vec_t> vecs_t;

std::map<std::pair<size_t, size_t>, vecs_t> tup_bank;

vecs_t &gen_tup(size_t s, size_t n) {
  auto it = tup_bank.find(std::make_pair(s, n));
  if (it != tup_bank.end()) {
    return it->second;
  }
  tup_bank[std::make_pair(s, n)] = vecs_t();
  auto &res = tup_bank[std::make_pair(s, n)];
  if (s == 0) {
    vec_t new_p(n);
    std::fill(new_p.begin(), new_p.end(), 0);
    res.push_back(new_p);
  } else if (n == 1) {
    res.push_back({s});
  } else {
    for (size_t i = 0; i <= s; i++) {
      for (auto p : gen_tup(s - i, n - 1)) {
        auto new_p(p);
        new_p.push_back(i);
        res.push_back(new_p);
      }
    }
  }
  return res;
}

typedef std::multiset<size_t> row_t;
typedef std::set<row_t> rows_t;

const size_t VEC_SIZE = 6;
const vec_t PRIMES{2, 3, 5, 7, 11, 13, 17, 19, 23};

rows_t gen_rows(vec_t &pows) {
  rows_t rows, new_rows;
  row_t init;
  for (size_t i = 0; i < VEC_SIZE; i++) {
    init.insert(1);
  }
  rows.insert(init);

  for (size_t i = 0; i < pows.size(); i++) {
    auto p = PRIMES[i];
    auto pow = pows[i];
    auto ptups = gen_tup(pows[i], VEC_SIZE);

    new_rows.clear();
    for (auto row : rows) {
      for (auto tup : ptups) {
        row_t new_row;
        auto r = tup.begin();
        for (auto elt : row) {
          auto new_elt = elt;
          for (size_t j = 0; j < *r; j++) {
            new_elt *= p;
          }
          new_row.insert(new_elt);
          r++;
        }
        new_rows.insert(new_row);
      }
    }
    std::swap(rows, new_rows);
  }

  return rows;
}

typedef std::map<size_t, std::vector<row_t>> rowdict_t;

rowdict_t rows_to_rowdict(rows_t all_rows) {
  rowdict_t rowdict;

  // split by sum
  auto it = all_rows.begin();
  while (it != all_rows.end()) {
    if (std::adjacent_find(it->begin(), it->end()) == it->end()) {
      size_t sum = 0;
      for (auto elt : *it) {
        sum += elt;
      }
      if (rowdict.find(sum) == rowdict.end()) {
        rowdict[sum] = {};
      }
      rowdict[sum].push_back(*it);
    }
    it = all_rows.erase(it);
  }

  size_t total = 0;
  for (auto pair : rowdict) {
    total += pair.second.size();
  }
  std::cout << ">>> (enum) total possible 6-vecs: " << total << "\n";

  // reduction
  auto rows_it = rowdict.begin();
  while (rows_it != rowdict.end()) {
    auto S = rows_it->first;
    auto rows = rows_it->second;

    // if < 12, remove and continue
    if (rows.size() < 12) {
      rows_it = rowdict.erase(rows_it);
      continue;
    }

    // counter
    std::map<size_t, size_t> count;
    for (auto row : rows) {
      for (auto elt : row) {
        if (count.find(elt) == count.end()) {
          count[elt] = 0;
        }
        count[elt]++;
      }
    }

    // remove while some count[elt] is 1
    bool did_remove = true;
    while (did_remove && rows.size() >= 12) {
      did_remove = false;
      auto it = rows.begin();
      while (it != rows.end()) {
        if (std::any_of(it->begin(), it->end(),
                        [&count](size_t elt) { return count[elt] == 1; })) {
          for (auto elt : *it) {
            count[elt]--;
          }
          did_remove = true;
          it = rows.erase(it);
        } else {
          it++;
        }
      }
    }

    // if < 12, remove and continue
    if (rows.size() < 12) {
      rows_it = rowdict.erase(rows_it);
      continue;
    } else {
      rows_it->second = rows;
      rows_it++;
    }
  }

  return rowdict;
}

int main(int argc, char *argv[]) {
  char *filename = argv[1];
  vec_t exponents(argc - 2);

  for (size_t i = 0; i < exponents.size(); i++) {
    exponents[i] = std::stol(argv[i + 2], NULL, 10);
  }

  auto rows = rows_to_rowdict(gen_rows(exponents));

  // write rows
  std::ofstream file;
  file.open(filename);

  uint64_t nrows = 0;

  for (auto vv : rows) {
    nrows += vv.second.size();
    for (auto v : vv.second) {
      auto it = v.rbegin();
      file << *it++;
      while (it != v.rend()) {
        file << " " << *it++;
      }
      file << "\n";
    }
  }

  std::cout << ">>> (enum) wrote " << nrows << " rows to " << filename
            << " in ascending order of sum\n";
}
