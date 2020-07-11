/*
 * SampleDecoder.cpp
 *
 * For more information, see SampleDecoder.h
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <cfloat>
#include <limits>

#include "VSBPPDecoder.h"
#include "global.h"

using namespace std;

#define TESTE 0

#define DEBUG 0
#define CHANGE_CHROM 1

SampleDecoder::SampleDecoder(char *instance) //constructor
{
    /* realiza a leitura dos dados de entrada e inicializa os atributos da classe */
    read_prob_sizes(instance);
    binCap = (int *) malloc(sizeof(int) * bin_types);
    binCost = (int *) malloc(sizeof(int) * bin_types);
    weight = (int *) malloc(sizeof(int) * items);
	bin_item = (int *) malloc(sizeof(int) * items+1); /* inclui no ficticio */
    read_instance(instance);
}

SampleDecoder::~SampleDecoder() {
    free(binCap);
    free(binCost);
    free(weight);
	free(bin_item);
    binCap  = NULL;
    binCost = NULL;
    weight  = NULL;
	bin_item= NULL;
} // destructor

double SampleDecoder::decode(BRKGA::Chromosome& chromosome, bool /* not-used */) const {

    double myFitness = 0.0;
    int k, z;
	char key[3000] = { NULL };
	char value[3000] = { NULL }; // key => chromossome, value => mapped int vector
	

    std::map<unsigned,double> myMap;
    typedef std::pair< double, unsigned > ValueKeyPair;
    std::vector< ValueKeyPair > rank(chromosome.size());
    
    
    for (unsigned i = 0; i < chromosome.size(); ++i) {
        rank[i] = ValueKeyPair(chromosome[i], i);
        myMap[i] = chromosome[i];
    }
    
	std::sort(rank.begin(), rank.end());
	vector<int> individuo;
	for (std::vector< ValueKeyPair >::const_iterator i = rank.begin(); i != rank.end(); ++i) {
		individuo.push_back(i->second);
	}

    myFitness = calc_fitness(individuo);
#if DEBUG
	for(k = 0; k < items; k++)
		printf("%.2f ", chromosome[k]);
	printf("%.0f\n", myFitness);
#endif
	
	for(k = 0; k < items; k++)
	{
		char tmp_k[50] = { NULL };
		char tmp_v[50] = { NULL };
		sprintf(tmp_k, "%.2f;", chromosome[k]);
		sprintf(tmp_v, "%d;", individuo[k]);
		if(k == 0)
		{
			strcpy(key, tmp_k);
			strcpy(value, tmp_v);
		} else {
			strcat(key, tmp_k);
			strcat(value, tmp_v);
		}
	}
	string kk(key), vv(value);
	//global::hSolutions[kk] = vv;
	if(global::hSolutions.size() <= 50 || (global::hSolutions.find(kk) == global::hSolutions.end() && myFitness > 1.05*global::bestsol) )
	{
		global::hSolutions[kk] = vv;
		myFitness = tabu_search(individuo, myFitness);
	} /*else {
		cout << "TS nao executado: " << myFitness << endl;
	}
	*/
#if CHANGE_CHROM
    for(k = 0; k < items; k++)
		chromosome[k] = myMap[individuo[k]];
#endif
	memset(key, 0, 3000);
	memset(value, 0, 3000);

	for(k = 0; k < items; k++)
	{
		char tmp_k[50] = { NULL };
		char tmp_v[50] = { NULL };
		sprintf(tmp_k, "%.2f;", chromosome[k]);
		sprintf(tmp_v, "%d;", individuo[k]);
		if(k == 0)
		{
			strcpy(key, tmp_k);
			strcpy(value, tmp_v);
		} else {
			strcat(key, tmp_k);
			strcat(value, tmp_v);
		}
	}

	string kk2(key), vv2(value);
	global::hSolutions[kk2] = vv2;

#if DEBUG
	myFitness = calc_fitness(individuo);
	for(k = 0; k < items; k++)
		printf("%.2f ", chromosome[k]);
	printf("%.0f\n", myFitness);

	for(k = 0; k < items; k++)
		printf("%d ", individuo[k]);
	printf("%.0f\n", myFitness);

	printf("\n");
#endif
	if(myFitness < global::bestsol) global::bestsol = myFitness;
	return myFitness;
}

int SampleDecoder::getnbObj(){ return items; }

/*****************************************************************************/
double SampleDecoder::calc_fitness(vector<int> & individuo) const
{
    double f = 0.0;
    int *Q; /* Peso acumulado */
    int *etiqueta_custo;
    int *etiqueta_pred;
    int i, j;
    int custo_acumulado = 0;
	int peso = 0;
	double penalidade = 0.0;
	double fi_bi = DBL_MAX; /* minima razao entre fi / bi (custo / capacidade) */
	double tmp = 0.0;
	int peso_excesso = 0;
	int peso_excesso_min = INT_MAX;
    int binNumber = 1;

    Q = (int *) malloc(sizeof(int) * (items + 1));
    etiqueta_custo = (int *) malloc(sizeof(int) * (items + 1));
    etiqueta_pred = (int *) malloc(sizeof(int) * (items + 1));
    
    for(i = 0; i < items + 1; i++)
    {
        etiqueta_custo[i] = 0;
        etiqueta_pred[i] = -1;
    }
    
    Q[0] = 0;
    for(i = 1; i < items + 1; i++)
        Q[i] = Q[i - 1] + weight[individuo[i - 1]];

    for(i = 0; i < items; i++)
    {
        for(j = i + 1; j < items+1; j++)
        {
            int val = Q[j] - Q[i];
            int entrou = 0;
            custo_acumulado = etiqueta_custo[i];
            if(val <= binCap[0])
            {
                custo_acumulado += binCost[0];
                entrou = 1;
            }
            else if(val <= binCap[1])
            {
                custo_acumulado += binCost[1];
                entrou = 1;
            }
            else if(val <= binCap[2])
            {
                custo_acumulado += binCost[2];
                entrou = 1;
            }
            
            if(entrou)
            {
                if(etiqueta_pred[j] == -1 || custo_acumulado < etiqueta_custo[j])
                {
                    etiqueta_custo[j] = custo_acumulado /*+ etiqueta_custo[i]*/;
					etiqueta_pred[j] = i;
				}
			} else break;
		}
	}

	f = etiqueta_custo[items];

	i = items;
	bin_item[i] = binNumber;
	do
	{
		for(j = i-1; j >= etiqueta_pred[i]; j--)
			if(j < items)
				bin_item[j] = binNumber;

		i = etiqueta_pred[i];
		binNumber++;
	} while(i > 0);

	free(Q);
	free(etiqueta_custo);
	free(etiqueta_pred);
	Q = NULL;
	etiqueta_custo = NULL;
	etiqueta_pred = NULL;
    return f;
}

/*****************************************************************************/
void SampleDecoder::read_prob_sizes(char * arquivo)
{
    FILE *fp;
    int tmp = 0;
    fp = fopen(arquivo, "r");
    if(fp == NULL)
        exit(1);
    
    fscanf(fp, "%d", &tmp);
    items = tmp;
    fscanf(fp, "%d", &tmp);
    bin_types = tmp;
    
    fclose(fp);
}

/*****************************************************************************/
void SampleDecoder::read_instance(char *arquivo)
{
    FILE *fp;
    int i = 0;
    int j = 0;
    int tmp = 0;
    int *d;
    
    d = (int *) malloc(sizeof(int) * bin_types*2);
    
    fp = fopen(arquivo, "r");
    if(fp == NULL)
        exit(1);
    
    for(i=0; i<4; i++)
        fscanf(fp, "%d", &tmp);
    
    for(j = 0; j < bin_types*2; j++)
    {
        fscanf(fp, "%d", &tmp);
        d[j] = tmp;
    }
    
    i = 0;
    for(j = 0; j < bin_types; j++)
    {
        binCap[j] = d[i];
        binCost[j] = d[i+1];
        i+=2;
    }
    
    for(i=0; i<items; i++)
    {
        fscanf(fp, "%d", &tmp);
        weight[i] = tmp;
    }
    
    free(d);
    fclose(fp);
}

/******************************************************************************/
/******************************************************************************/
void SampleDecoder::best_move(move *mov, int iter, double c_best, double c_curr,
               int **tabu_time, vector<int> & individuo) const
{
	int i, j;
	mov->value = DBL_MAX;

	for(i = 0; i < items; i++)
	{
		for(j = 0; j < items; j++)
		{
			if(bin_item[i] != bin_item[j] && weight[individuo[i]] != weight[individuo[j]])
			{
				int tmp = individuo[i];
				individuo[i] = individuo[j];
				individuo[j] = tmp;

				c_curr = calc_fitness(individuo);
				if(tabu_time[i][j] < iter || c_curr < c_best)
				{
					if(c_curr < mov->value)
					{
						mov->value = c_curr;
						mov->i = i;
						mov->j = j;
					}
				}
				tmp = individuo[j];
				individuo[j] = individuo[i];
				individuo[i] = tmp;
			}
		}
	}
}
/******************************************************************************/
void SampleDecoder::execute_move(move *mov, vector<int> & individuo) const
{
    int tmp = individuo[mov->i];
    individuo[mov->i] = individuo[mov->j];
    individuo[mov->j] = tmp;
}
/******************************************************************************/
double SampleDecoder::tabu_search(vector<int> & individuo, double cost) const
{
    double ts_cost = 0.0;
    int i, j;
    int iter, max_iter, tabu_size;
    int **tabu_time;
    vector<int> tmp_sol(items+1);
    int *best_sol;
    double c_best, c_curr;
    move mov;
    tabu_time = imatrix(0, items, 0, items);
    best_sol = (int *) malloc(sizeof(int) * items+1);
    
    
    for(i = 0; i < items; i++)
    {
        tmp_sol[i] = individuo[i];
        best_sol[i] = individuo[i];
        for(j = 0; j < items; j++)
            tabu_time[i][j] = 0;
    }
    tmp_sol[items] = individuo[items];
    best_sol[items]= individuo[items];
    
    iter = 0;
    max_iter = 10; /* Quantidade de iteracoes a serem executadas no TS */
    tabu_size = 5;  /* Tabu tenure */
    c_curr = cost;
    c_best = c_curr;
    
    while(iter < max_iter)
    {
        iter++;
        best_move(&mov, iter, c_best, c_curr, tabu_time, tmp_sol);
		execute_move(&mov, tmp_sol);
        
        tabu_time[mov.i][mov.j] = iter + tabu_size;
        
        c_curr = mov.value;
        if(c_curr < c_best)
        {
            c_best = c_curr;
			individuo = tmp_sol;
        }
    }
    ts_cost = c_best;

    free_imatrix(tabu_time, 0, items, 0);
    free(best_sol);
    return ts_cost;
}
/******************************************************************************/
int ** SampleDecoder::imatrix(int nrl, int nrh, int ncl, int nch) const
{
    int i;
    int **m;
    m = (int **)malloc((unsigned)(nrh - nrl +1)*sizeof(int*));
    if (!m)
    {
        fprintf(stderr,"\n%s\n", "allocation failure 1 in imatrix()");
        fprintf(stderr,"... now exiting to system ...\n");
        exit(1);
    }
    m -= nrl;
    for(i = nrl; i <= nrh ; i++)
    {
        m[i]= (int*)malloc((unsigned)(nch - ncl +1)*sizeof(int));
        if (!m[i])
        {
            fprintf(stderr,"\n%s\n", "allocation failure 2 in imatrix()");
            fprintf(stderr,"... now exiting to system ...\n");
            exit(1);
        }
        m[i] -= ncl;
    }
    return m;
}
/******************************************************************************/
void SampleDecoder::free_imatrix(int **m, int nrl, int nrh, int ncl) const
{
    int i;
    for (i = nrh; i >= nrl; i--)
        free((char*) (m[i] + ncl));
    free((char *)(m + nrl));
}
/******************************************************************************/

