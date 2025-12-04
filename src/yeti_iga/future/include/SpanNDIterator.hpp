#pragma once
#include "BSplineTensor.hpp"

#include <iostream>  // Temp for debug

class SpanNDIterator {
public:
    SpanNDIterator(const BSplineTensor& tensor_)
        : tensor(tensor_), n_dims(tensor.components.size()), done(false)
    {
        // Build list of valid spans for each dimension
        valid_spans.resize(n_dims);
        std::cout << "ndims : " << n_dims << "\n";
        for (size_t d = 0; d < n_dims; ++d) {
            std::cout << "d : " << d << "\n";
            const auto& kv = tensor.components[d].getKnotVector();
            int deg = tensor.components[d].getDegree();
            int n = kv.size() - 1;
            std::cout << "n : " << n << "\n";


            for (int i = deg; i <= n - deg - 1; ++i) {
                std::cout << "boucle for i\n";
                if (kv[i+1] - kv[i] > 0.)
                    valid_spans[d].push_back(i);
            }
        }

        // Initialize current index (zero for all)
        idx.resize(n_dims, 0);
        done = (n_dims == 0 || valid_spans[0].empty());
        std::cout << valid_spans[0].size() << "\t" << valid_spans[1].size() << "\n";
        std::cout << done << "\n";
    }

    const std::vector<int> current() const {
        std::vector<int> spans(n_dims);
        for (size_t d = 0; d < n_dims; ++d)
            spans[d] = valid_spans[d][idx[d]];
        return spans;
    }

    void next() {
        if (done) return;
        for (int d = n_dims-1; d >= 0; --d) {
            if (++idx[d] < valid_spans[d].size())
                return;
            idx[d] = 0;
        }
        done = true;
    }

    bool is_done() const {return done;}

private:
    const BSplineTensor& tensor;
    size_t n_dims;
    std::vector<std::vector<int>> valid_spans;
    std::vector<size_t> idx;
    bool done;
};