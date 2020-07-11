/**
 * SampleDecoder.h
 *
 * The chromosome's fitness is computed as the sum of its individual alleles. We also show how to
 * generate a random permutation of {0, 1, ..., n-1} using the supplied chromosome.
 *
 * Any Decoder class must implement either of the methods:
 *     1. double decode(std::vector< double >&)
 *     2. double decode(std::vector< double >&) const
 *     3. double decode(const std::vector< double >&)
 *     4. double decode(const std::vector< double >&) const
 * where the returned double corresponds to the fitness of that chromosome. If parallel decoding
 * is to be used in the BRKGA framework, then decode() *must* be thread-safe; the best way to
 * guarantee this is by via member const-correctness -- see (2) and (4) above --  so that the
 * property will be checked at compile time. An exception to this rule is the use of the mutant
 * modifier in data members of the decoder, which bypasses const-correctness and may introduce
 * dangerous race conditions in a multithreaded environment.
 *
 * The chromosome inside the BRKGA framework can be changed if desired. To do so, just use the
 * signatures (1) and (3) above which do not forbid the chromosome to be changed (and such change
 * gets reflected in BRKGA). Please use double values in the interval [0,1) when updating, thus
 * obeying the BRKGA guidelines -- for speed, BRKGA does *not* enforce this.
 *
 * Created on : Nov 17, 2011 by rtoso
 * Authors    : Rodrigo Franco Toso <rtoso@cs.rutgers.edu>
 *              Mauricio G.C. Resende <mgcr@research.att.com>
 * Copyright 2010, 2011 Rodrigo Franco Toso and Mauricio G.C. Resende.
 *
 * This file is part of the BRKGA API.
 *
 * The BRKGA API is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The BRKGA API is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with the BRKGA API. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SAMPLEDECODER_H
#define SAMPLEDECODER_H

#include <list>
#include <vector>
#include <algorithm>
#include <string>

#include <stdio.h>
#include <limits.h>
#include "chromosome.hpp"

class SampleDecoder {
public:
    typedef struct move_info
    {
        double value;
        int i;
        int j;
    } move;
    
    SampleDecoder(char *instance);	// Constructor
    ~SampleDecoder();	// Destructor
    
    // Decode a chromosome, returning its fitness as a double-precision floating point:
    double decode(BRKGA::Chromosome& chromosome, bool rewrite) const;
    int getnbObj();
    double calc_fitness(std::vector<int> & individuo) const;
    void read_prob_sizes(char * arquivo);
    void read_instance(char *arquivo);

	void best_move(move *mov, int iter, double c_best, double c_curr, int **tabu_time, std::vector<int> & individuo) const;
    void execute_move(move *mov, std::vector<int> & individuo) const;
    double tabu_search(std::vector<int> & individuo, double cost) const;
    void nrerror(char *error_text);
    void free_imatrix(int **m, int nrl, int nrh, int ncl) const;
    int **imatrix(int nrl, int nrh, int ncl, int nch) const;
    
private:
    int bin_types;
    int items;
    int * binCap;
    int * binCost;
    int * weight;
	int * bin_item; /* numero do bin em que o item foi inserido */
};

#endif
