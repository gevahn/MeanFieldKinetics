#include <vector>
#include <math.h>
#include <string>
#include <fstream>
#include <iostream>
#include <random>
#include <algorithm>
#include <iterator>
#include <functional>
#include <assert.h>
#include "tclap/CmdLine.h"

using namespace TCLAP;
using namespace std;

#define PI 3.14159265
#define h_bar 6.58211928e-16	//eV.s 
#define K_B 8.6173324e-5		//eV.K-1
#define T 300					//K



class KMC {

public:

	enum Species { Hatom, H2mol, Vacancy };

	struct Site {
		void (KMC::*reaction)(int, int);
		Species type;
		double rate;
	};


	vector< vector < Site* >> getNN(int i, int j){
		vector<vector<Site*>> NNlist(3);
		
		int iplus = (i + 1) % numberOfSites;
		int iminus = (i - 1 + numberOfSites) % numberOfSites;
		int jplus = (j + 1) % numberOfSites;
		int jminus = (j - 1 + numberOfSites) % numberOfSites;

		NNlist[grid[iplus][j].type].push_back(&grid[iplus][j]);
		NNlist[grid[iplus][jplus].type].push_back(&grid[iplus][jplus]);
		NNlist[grid[iplus][jminus].type].push_back(&grid[iplus][jminus]);
		NNlist[grid[i][jplus].type].push_back(&grid[i][jplus]);
		NNlist[grid[i][jminus].type].push_back(&grid[i][jminus]);
		NNlist[grid[iminus][j].type].push_back(&grid[iminus][j]);
		NNlist[grid[iminus][jplus].type].push_back(&grid[iminus][jplus]);
		NNlist[grid[iminus][jminus].type].push_back(&grid[iminus][jminus]);

		return NNlist;
	}

	void Hreaction(int i, int j) {
		double rateSoFar = kd;
		double randRate;
		randRate= dist(mersenneEngine)*totalRate;

		// desorption
		if (randRate< rateSoFar) {
			grid[i][j].reaction = &KMC::Vacreaction;
			grid[i][j].type = Vacancy;
			numberOfH -= 1;
			numberOfVac += 1;
			return;
		}

		vector<vector<Site*>> NNlist;
		NNlist = getNN(i,j);		
		

		// diffusion
		rateSoFar += kdiff*NNlist[Vacancy].size();
		if (randRate< rateSoFar) {
			grid[i][j].reaction = &KMC::Vacreaction;
			grid[i][j].type = Vacancy;
			if (NNlist[Vacancy].size() - 1 < 0) cout << "meow1" <<endl;			
			uniform_int_distribution<int> dn(0,NNlist[Vacancy].size() - 1);	
			int chooseSite = dn(mersenneEngine);
			NNlist[Vacancy][chooseSite]->reaction = &KMC::Hreaction;
			NNlist[Vacancy][chooseSite]->type = Hatom;
			return;

		}

		// recombination
		rateSoFar += grid[i][j].rate * NNlist[Hatom].size();
		if (randRate< rateSoFar) {			
			grid[i][j].reaction = &KMC::Vacreaction;
			grid[i][j].type = Vacancy;
			uniform_int_distribution<int> dn(0,NNlist[Hatom].size() - 1);	
			int chooseSite = dn(mersenneEngine);
			NNlist[Hatom][chooseSite]->reaction = &KMC::H2reaction;
			NNlist[Hatom][chooseSite]->type = H2mol;

			numberOfH -= 2;
			numberOfH2 += 1;
			numberOfVac += 1;
			return;

		}
		numberOfNullSteps += 1;

		return;

	}		

	void H2reaction(int i, int j) {
                double rateSoFar = kr;
                double randRate;
                randRate= dist(mersenneEngine)*totalRate;

                // desorption
                if (randRate< rateSoFar) {
                        grid[i][j].reaction = &KMC::Vacreaction;
                        grid[i][j].type = Vacancy;
			numberOfVac += 1;
			numberOfH2 -= 1;
                        return;
                }

                vector<vector<Site*>> NNlist;
                NNlist = getNN(i,j);

/*
                // diffusion
                rateSoFar += kdiff*NNlist[Vacancy].size();
                if (randRate< rateSoFar) {
                        grid[i][j].reaction = Vac:raction;
                        grid[i][j].type = Vacancy;
                        uniform_int_distribution<int> dn(0,NNlist[Vacancy].size() - 1);
                        int chooseSite = dn(mersenneEngine);
                        NNlist[Vacancy][chooseSite]->reaction = H:reaction;
                        NNlist[Vacancy][chooseSite]->type = Hatom;

                }
*/
                // dissociation
                rateSoFar += grid[i][j].rate/Keq*NNlist[Vacancy].size();
                if (randRate< rateSoFar) {
                        grid[i][j].reaction = &KMC::Hreaction;
                        grid[i][j].type = Hatom;
                        uniform_int_distribution<int> dn(0,NNlist[Vacancy].size() - 1);
                        int chooseSite = dn(mersenneEngine);
                        NNlist[Vacancy][chooseSite]->reaction = &KMC::Hreaction;
                        NNlist[Vacancy][chooseSite]->type = Hatom;

			numberOfH2 -= 1;
			numberOfH += 2;
			numberOfVac -= 1;
			return;

                }
		numberOfNullSteps += 1;

                return;

        }


	void Vacreaction(int i, int j) {
                double rateSoFar = ka;
                double randRate;
                randRate= dist(mersenneEngine)*totalRate;

                // adsortion
                if (randRate< rateSoFar) {
                        grid[i][j].reaction = &KMC::Hreaction;
                        grid[i][j].type = Hatom;
			numberOfVac -= 1;
			numberOfH += 1;
                        return;
                }

                vector<vector<Site*>> NNlist;
                NNlist = getNN(i,j);


                // diffusion
                rateSoFar += kdiff*NNlist[Hatom].size();
                if (randRate< rateSoFar) {
                        grid[i][j].reaction = &KMC::Hreaction;
                        grid[i][j].type = Hatom;
                        uniform_int_distribution<int> dn(0,NNlist[Hatom].size() - 1);
                        int chooseSite = dn(mersenneEngine);
                        NNlist[Hatom][chooseSite]->reaction = &KMC::Vacreaction;
                        NNlist[Hatom][chooseSite]->type = Vacancy;

			return;
                }

                // dissociation
                rateSoFar += grid[i][j].rate/Keq*NNlist[H2mol].size();
                if (randRate< rateSoFar) {
                        grid[i][j].reaction = &KMC::Hreaction;
                        grid[i][j].type = Hatom;
                        uniform_int_distribution<int> dn(0,NNlist[H2mol].size() - 1);
                        int chooseSite = dn(mersenneEngine);
                        NNlist[H2mol][chooseSite]->reaction = &KMC::Hreaction;
                        NNlist[H2mol][chooseSite]->type = Hatom;

			numberOfVac -= 1;
			numberOfH += 2;
			numberOfH2 -= 1;
			return;
                }
		numberOfNullSteps += 1;
                return;

        }








	KMC(int _numberOfSites, double _ka, double _kd, double _kdiff, double _kr, double _Keq, double _mean, double _sigma, int _numberOfH, int _numberOfVac, int _numberOfH2) {
		numberOfSites = _numberOfSites;

		grid = vector< vector< Site>>(numberOfSites, vector<Site>(numberOfSites));
		ka = _ka;	
		kd = _kd;	
		kdiff = _kdiff;	
		kr = _kr;	
		Keq = _Keq;	
		mean = _mean;	
		sigma = _sigma;	

//		double rates[] = {8*kdiff+ka,kd+8*kdiff+8*kr
		totalRate = 8* (kdiff+ka+kd+kr+exp(mean));

		numberOfH = _numberOfH;	
		numberOfVac = _numberOfVac;	
		numberOfH2 = _numberOfH2;	
		numberOfNullSteps = 0;

		random_device rd;
		mersenneEngine = mt19937(rd());
		dist = uniform_real_distribution<double> (0.0,1.0);
		normal_distribution<double> nd(mean, sigma);
		
		int currentGridPoint = 0;	
	
		for (int assigned = 0; assigned < numberOfH; assigned++) {
			int i = currentGridPoint / numberOfSites;
			int j = currentGridPoint % numberOfSites;
			grid[i][j].reaction = &KMC::Hreaction;	 
			grid[i][j].rate =  exp(nd(mersenneEngine));	 
			grid[i][j].type = Hatom;	 
			currentGridPoint += 1;
		}
		for (int assigned = 0; assigned < numberOfH2; assigned++) {
			int i = currentGridPoint / numberOfSites;
			int j = currentGridPoint % numberOfSites;
			grid[i][j].reaction = &KMC::H2reaction;	 
			grid[i][j].rate =  exp(nd(mersenneEngine));	 
			grid[i][j].type = H2mol;	 
			currentGridPoint += 1;
		}
		for (int assigned = 0; assigned < numberOfVac; assigned++) {
			int i = currentGridPoint / numberOfSites;
			int j = currentGridPoint % numberOfSites;
			grid[i][j].reaction = &KMC::Vacreaction;	 
			grid[i][j].rate =  exp(nd(mersenneEngine));	 
			grid[i][j].type = Vacancy;	 
			currentGridPoint += 1;
		}
		for (int i =0; i < numberOfSites; i++){
			shuffle(begin(grid[i]), end(grid[i]), mersenneEngine);
		}
		shuffle(begin(grid), end(grid), mersenneEngine);

	}
	//rates
	double ka;
	double kd;
	double kdiff;
	double kr;
	double Keq;
	double mean;
	double sigma;

	double totalRate;

	int A;	
	int B;	
	int C;

	int numberOfNullSteps;

	int numberOfSites;	
        int numberOfH;
        int numberOfVac;
        int numberOfH2;
	
	mt19937 mersenneEngine;	
	uniform_real_distribution<double> dist;

	vector< vector <Site>> grid;

	void react(int i, int j) {
		(this->*grid[i][j].reaction)(i,j);
	}

	void KMCstep() {
		int i;
		int j;
		uniform_int_distribution<int> dn(0,numberOfSites - 1);
		i = dn(mersenneEngine);	
		j = dn(mersenneEngine);	
		react(i,j);
	}

	void print() {
		int n = numberOfSites * numberOfSites;
		cout << "Number of H: " << numberOfH << endl;
		cout << "Number of H2: " << numberOfH2 << endl;
		cout << "Number of Vac: " << numberOfVac << endl;
		cout << "\% of H: " << numberOfH * 1.0 / n << endl;
		cout << "\% of H2: " << numberOfH2 * 1.0 / n<< endl;
		cout << "\% of Vac: " << numberOfVac * 1.0 / n<< endl;
		cout << "Number of Null steps: " << numberOfNullSteps << endl;
	}

	void printGrid() {
		for(int i = 0; i < numberOfSites; i++) {
			for(int j = 0; j < numberOfSites; j++) {
				cout << grid[i][j].type << " ";
			}
			cout << endl;
		}
	}

	void resetGrid() {
		normal_distribution<double> nd(mean, sigma);
		

		for(int i; i < numberOfSites; i++) {
			for(int j; j < numberOfSites; j++) {
				grid[i][j].rate = nd(mersenneEngine);
				
			}
		}

	}



};



int main(int argc, char** argv) {

	int numberOfSites = 40;
	double ka;
	double kd = 0.1;	
	double kdiff = 1;
	double kr = 0;
	double Keq = 50;
	double mean = -6;
	double sigma = 0;
	
	try {
		CmdLine cmd("Enter rates and stuff");
		
		ValueArg<int> numberOfSitesArg("n", "numberOfSites", "adsorption rate", true,40,"int");
		cmd.add(numberOfSitesArg);
		ValueArg<double> kaArg("a", "ka", "adsorption rate", true,0.1,"double");
		cmd.add(kaArg);
		ValueArg<double> kdArg("d", "kd", "desorption rate", true,0.1,"double");
		cmd.add(kdArg);
		ValueArg<double> kdiffArg("i", "kdiff", "diffusion rate", true,1,"double");
		cmd.add(kdiffArg);
		ValueArg<double> krArg("r", "kr", "H2 desorption rate", true,0.1,"double");
		cmd.add(krArg);
		ValueArg<double> KeqArg("e", "Keq", "H2 equilibrium constant", true,50,"double");
		cmd.add(KeqArg);
		ValueArg<double> meanArg("m", "mean", "H2 reaction rate mean", true, -6 ,"double");
		cmd.add(meanArg);
		ValueArg<double> sigmaArg("s", "sigma", "H2 reaction rate sigma", true,0,"double");
		cmd.add(sigmaArg);

		cmd.parse(argc, argv);

		ka = kaArg.getValue();
		kd = kdArg.getValue();
		kdiff = kdiffArg.getValue();
		kr = krArg.getValue();
		Keq = KeqArg.getValue();
		mean = meanArg.getValue();
		sigma = sigmaArg.getValue();
		numberOfSites = numberOfSitesArg.getValue();
		

	} catch (ArgException &e) {
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		return 1;
	}


//	KMC(int _numberOfSites, double _ka, double _kd, double _kdiff, double _kr, double _Keq, double _Keq, double _sigma, int _numberOfH, int _numberOfVac, int _numberOfH2) {


	KMC sim = KMC(numberOfSites, ka, kd, kdiff, kr, Keq, mean, sigma, 0, numberOfSites * numberOfSites, 0);
	int numberOfMCSteps = 100000000;

	sim.print();
	cout << endl;
//	sim.printGrid();
	for (int i = 0; i < numberOfMCSteps; i++) {
		sim.KMCstep();
	}

	sim.print();
//	sim.printGrid();


	
	return 0;
}

