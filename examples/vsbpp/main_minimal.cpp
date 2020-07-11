/******************************************************************************
 * main_minimal.cpp: minimal code for calling BRKGA algorithms to solve
 *                   instances of the Traveling Salesman Problem.
 *
 * (c) Copyright 2015-2019, Carlos Eduardo de Andrade.
 * All Rights Reserved.
 *
 *  Created on : Mar 05, 2019 by andrade
 *  Last update: Mar 05, 2019 by andrade
 *
 * This code is released under LICENSE.md.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *****************************************************************************/


#include "VSBPPDecoder.h"
#include "brkga_mp_ipr.hpp"

#include <iostream>
#include <stdexcept>
#include <string>
#include <cstdio>
#include <limits>
using namespace std;

//-------------------------------[ Main ]------------------------------------//

int main(int argc, char* argv[]) {
    double tempo;
    clock_t start = clock();

    if(argc < 4) {
        cerr << "Usage: "<< argv[0]
             << " <seed> <config-file> <num-generations>"
                " <vsbpp-instance-file>" << endl;
        return 1;
    }

    try {
        ////////////////////////////////////////
        // Read command-line arguments and the instance
        ////////////////////////////////////////
        char instance_file[80];

        const unsigned seed = stoi(argv[1]);
        const string config_file = argv[2];
        const unsigned num_generations = stoi(argv[3]);
        sscanf(argv[4], "%s", instance_file);

        ////////////////////////////////////////
        // Read algorithm parameters
        ////////////////////////////////////////

        cout << "Reading parameters..." << endl;

        // C++14 syntax.
        auto params = BRKGA::readConfiguration(config_file);
        auto& brkga_params = params.first;

        // C++17 syntax.
        // auto [brkga_params, control_params] =
        //     BRKGA::readConfiguration(config_file);

        ////////////////////////////////////////
        // Build the BRKGA data structures and initialize
        ////////////////////////////////////////

        cout << "Building BRKGA data and initializing..." << endl;

        VSBPPDecoder decoder(instance_file);				// initializes the decoder and reads data
        const unsigned n = decoder.getnbObj();		// size of chromosomes = Number of alleles per chromosome

        BRKGA::BRKGA_MP_IPR<VSBPPDecoder> algorithm(
                decoder, BRKGA::Sense::MINIMIZE, seed,
                n, brkga_params);

        // NOTE: don't forget to initialize the algorithm.
        algorithm.initialize();

        ////////////////////////////////////////
        // Find good solutions / evolve
        ////////////////////////////////////////

        cout << "Evolving " << num_generations << " generations..." << endl;
        algorithm.evolve(num_generations);
	tempo = (double)(clock() - start)/CLOCKS_PER_SEC;

        auto best_cost = algorithm.getBestFitness();
        cout << "Best cost: " << best_cost << endl;
	printf("Time: %.4f\n", tempo);

    }
    catch(exception& e) {
        cerr << "\n***********************************************************"
             << "\n****  Exception Occured: " << e.what()
             << "\n***********************************************************"
             << endl;
        return 70; // BSD software internal error code
    }
    return 0;
}
