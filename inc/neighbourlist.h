#ifndef _NEIGHBOURLIST_H_
#define _NEIGHBOURLIST_H_

    #include <vector>

    namespace neigh {

        // Globals
        extern std::vector<unsigned int> adjncy;
        extern std::vector<double> Jijy_prime;
        extern std::vector<double> Jijz_prime;
        extern std::vector<double> Jijx_prime;
        extern std::vector<unsigned int> x_adj;
        extern std::vector<int> simspin;
        extern int nsimspin;

        // Functions
        void init();
        void ReadFile();
        void InteractionMatrix();
        void Heun(double Thermal_Fluct);

    }

#endif
