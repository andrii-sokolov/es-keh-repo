/*************************************/
/*
/*************************************/

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <omp.h>

/******************/
/*Loading elements of BOOST Library
/*(To have adaptive ODE solver)
/*NOTE: To compilt this code you need command:
/* g++ -O3 -o  oscillator_exe main_prime.cc -lboost_system -lboost_thread
/******************/
#include <boost/array.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

/******************/
/*GLOBAL VARIABLES*/
/******************/

//The mass of resonator [kg]
double m = 59e-6;
//Linear spring coefficient [N/m]
double k = 18.0;
//3-power [N/m^3]
double k3 = 0.0;
//Stopper distance [m]
double dstop = 47.5e-6;
//Stopper smooth delta [m]
double delta_dstop = 2e-6;
//Stopper elasticity
double kstop = 1.0e3;
//3-power spring coefficient [N/m^3]
double k3stop = 0.0;
//Stopper quality factor
double Qstop = 0.7;
//Stopper damping coefficient
double cstop = sqrt(m*kstop)/Qstop;
//Q-factor
double Q = 5.5;
//damping coefficient
double c = sqrt(m*k)/Q;
//External frequency
double fext = 10.0;
//External circular frequancy
double wext = 2*M_PI*fext;
//Amplitude
double Aext = 9.81*0.5;

//Initial C0 [F]
double C0 = 3.2e-12;
double Cpar = 28.0e-12;
//Load Resistance [Ohm]
double Rl = 6.65e6;
//Bias Voltage [V]
double V0 = 26.0;
//Capacitor distance [m]
double dcap = 54.2e-6;
// Numerical delta_t
double solv_delta = 0.01;

//Type to hold the data of the ODE solver
typedef boost::array< double , 3 > state_type;

// The type of the stepper parameters
typedef controlled_runge_kutta< runge_kutta_cash_karp54< state_type > > stepper_type;

// The relative error of ODE solver
double rel_err = 1.0e-17;
stepper_type controlled_stepper( rel_err = 1.0e-17);

//kstop
double kstop_function(double x)
{
	double k = 0.0;
	double x1 = dstop - delta_dstop;
	double x2 = dstop + delta_dstop;
	if(abs(x)<x1)
		k = 0.0;
	if((x>=x1)&(x<x2))
		k = ((kstop*pow(x1,3)-3.0*kstop*pow(x1,2)*x2)/(-pow(x2,3) + 3.0*x1*x2*x2 - 3.0*x1*x1*x2 + pow(x1,3))) +
				(6.0*kstop*x1*x2/(-pow(x2,3)+3.0*x1*x2*x2-3.0*x1*x1*x2 + pow(x1,3)))*x +
				-((3.0*kstop*x2 + 3.0*kstop*x1)/(-pow(x2,3)+3.0*x1*x2*x2-3.0*x1*x1*x2 + pow(x1,3)))*x*x +
				(2.0*kstop/(-pow(x2,3)+3.0*x1*x2*x2 -3.0*x1*x1*x2 + pow(x1,3)))*x*x*x;
	if(x>=x2)
		k = kstop;
	if((x<=-x1)&(x>-x2))
		k = ((kstop*pow(x1,3)-3.0*kstop*pow(x1,2)*x2)/(-pow(x2,3) + 3.0*x1*x2*x2 - 3.0*x1*x1*x2 + pow(x1,3))) +
				-(6.0*kstop*x1*x2/(-pow(x2,3)+3.0*x1*x2*x2-3.0*x1*x1*x2 + pow(x1,3)))*x +
				-((3.0*kstop*x2 + 3.0*kstop*x1)/(-pow(x2,3)+3.0*x1*x2*x2-3.0*x1*x1*x2 + pow(x1,3)))*x*x +
				-(2.0*kstop/(-pow(x2,3)+3.0*x1*x2*x2 -3.0*x1*x1*x2 + pow(x1,3)))*x*x*x;
	if(x<=-x2)
		k = kstop;
	return(k);
}

double k3stop_function(double x)
{
	double k = 0.0;
	double x1 = dstop - delta_dstop;
	double x2 = dstop + delta_dstop;
	if(abs(x)<x1)
		k = 0.0;
	if((x>=x1)&(x<x2))
		k = ((k3stop*pow(x1,3)-3.0*k3stop*pow(x1,2)*x2)/(-pow(x2,3) + 3.0*x1*x2*x2 - 3.0*x1*x1*x2 + pow(x1,3))) +
				(6.0*k3stop*x1*x2/(-pow(x2,3)+3.0*x1*x2*x2-3.0*x1*x1*x2 + pow(x1,3)))*x +
				-((3.0*k3stop*x2 + 3.0*k3stop*x1)/(-pow(x2,3)+3.0*x1*x2*x2-3.0*x1*x1*x2 + pow(x1,3)))*x*x +
				(2.0*k3stop/(-pow(x2,3)+3.0*x1*x2*x2 -3.0*x1*x1*x2 + pow(x1,3)))*x*x*x;
	if(x>=x2)
		k = k3stop;
	if((x<=-x1)&(x>-x2))
		k = ((k3stop*pow(x1,3)-3.0*k3stop*pow(x1,2)*x2)/(-pow(x2,3) + 3.0*x1*x2*x2 - 3.0*x1*x1*x2 + pow(x1,3))) +
				-(6.0*k3stop*x1*x2/(-pow(x2,3)+3.0*x1*x2*x2-3.0*x1*x1*x2 + pow(x1,3)))*x +
				-((3.0*k3stop*x2 + 3.0*k3stop*x1)/(-pow(x2,3)+3.0*x1*x2*x2-3.0*x1*x1*x2 + pow(x1,3)))*x*x +
				-(2.0*k3stop/(-pow(x2,3)+3.0*x1*x2*x2 -3.0*x1*x1*x2 + pow(x1,3)))*x*x*x;
	if(x<=-x2)
		k = k3stop;
	return(k);
}

//cstop
double cstop_function(double x)
{
	double c = 0.0;
	double x1 = dstop - delta_dstop;
	double x2 = dstop + delta_dstop;
	if(abs(x)<x1)
		c = 0.0;
	if((x>=x1)&(x<x2))
		c = ((cstop*pow(x1,3)-3.0*cstop*pow(x1,2)*x2)/(-pow(x2,3) + 3.0*x1*x2*x2 - 3.0*x1*x1*x2 + pow(x1,3))) +
				(6.0*cstop*x1*x2/(-pow(x2,3)+3.0*x1*x2*x2-3.0*x1*x1*x2 + pow(x1,3)))*x +
				-((3.0*cstop*x2 + 3.0*cstop*x1)/(-pow(x2,3)+3.0*x1*x2*x2-3.0*x1*x1*x2 + pow(x1,3)))*x*x +
				(2.0*cstop/(-pow(x2,3)+3.0*x1*x2*x2 -3.0*x1*x1*x2 + pow(x1,3)))*x*x*x;
	if(x>=x2)
		c = cstop;
	if((x<=-x1)&(x>-x2))
		c = ((cstop*pow(x1,3)-3.0*cstop*pow(x1,2)*x2)/(-pow(x2,3) + 3.0*x1*x2*x2 - 3.0*x1*x1*x2 + pow(x1,3))) +
				-(6.0*cstop*x1*x2/(-pow(x2,3)+3.0*x1*x2*x2-3.0*x1*x1*x2 + pow(x1,3)))*x +
				-((3.0*cstop*x2 + 3.0*cstop*x1)/(-pow(x2,3)+3.0*x1*x2*x2-3.0*x1*x1*x2 + pow(x1,3)))*x*x +
				-(2.0*cstop/(-pow(x2,3)+3.0*x1*x2*x2 -3.0*x1*x1*x2 + pow(x1,3)))*x*x*x;
	if(x<=-x2)
		c = cstop;
	return(c);
}

//The stopper force with inelastic collision
double Fstop(double x, double v)
{
	double f = -kstop_function(x)*x - k3stop_function(x)*pow(x,3) - cstop_function(x)*v;
	return(f);
}

//Variable capacitor formula
double Cap(double x)
{
	return(Cpar + 2.0*pow(dcap,2)*C0/(pow(dcap,2) - pow(x,2)));
}

//Transducer force
double Ftrans(double x, double q)
{
	return(2.0*pow(q,2)*C0*pow(dcap,2)*x/(pow(pow(dcap,2)-pow(x,2),2)*pow(Cap(x),2)));
}

//Oscillation system of equations
void oscillator( const state_type &x , state_type &dxdt , double t )
{
    dxdt[0] = -k*x[1]/m - Ftrans(x[1],x[2])/m + Fstop(x[1],x[0])/m - c*x[0]/m + Aext*cos(wext*t);
    dxdt[1] = x[0];
    dxdt[2] = V0/Rl - x[2]/(Cap(x[1])*Rl);
}

//Structure to send data from ODE
struct push_back_state_and_time
{
    std::vector< state_type >& m_states;
    std::vector< double >& m_times;

    push_back_state_and_time( std::vector< state_type > &states , std::vector< double > &times )
    : m_states( states ) , m_times( times )
    {
    	m_times = {};
        m_states = {};
    }

    void operator()( const state_type &x , double t )
    {
        m_states.push_back( x );
        m_times.push_back( t );
    }
};

//Reader of the experimental data
void read_exp_data(std::string filename, vector<double> &freq, vector<double> &W)
{
	freq = {};
	W = {};
  	std::string line;
  	ifstream myfile (filename);

  	vector<string> temp;

  	if (myfile.is_open())
  	{
    	while ( getline (myfile,line) )
    	{
    		boost::split(temp, line, boost::is_any_of(","));
    		freq.push_back( stod(temp[0]));
    		W.push_back(stod(temp[1]));
    	}
    myfile.close();
  	}

  else cout << "Unable to open file";
}

//differennce between an experiment and a theoretical data
double dif(vector<double> W1, vector<double> W2)
{
	double sumdif = 0.0;
	for(int i=0; i<W1.size(); ++i)
		sumdif += abs(W2[i]-W1[i]);
	return(sumdif);
}

double compare( std::string expfilename )
{
   	// vectors to hold the experimental data
	vector<double> exp_f;
	vector<double> exp_W;
	read_exp_data(expfilename, exp_f, exp_W);

   	// vector to hold the theoretical solution
   	vector<double> th_W={};

   	// initial conditions
    state_type x0 = { 0.0 , 0.0, C0*V0};

    vector<state_type> x_vec;
	vector<double> times;

	double tmax;
	double dt;
	int npoints;

	double xmax;
	double xmin;

	double W;
	int i;
	double rolling_t;

   	for(int j=0; j<exp_f.size(); ++j)
    {
    	// external circular frequency
    	wext = 2*M_PI*exp_f[j];

    	// maximal time
    	tmax = 200.0*2*M_PI/wext;
    	dt = solv_delta*2*M_PI/wext;

	    // RK4 method
	    npoints = integrate_adaptive(controlled_stepper,oscillator,
	    	x0 , 0.0 , tmax , dt, push_back_state_and_time(x_vec, times));

	    x0[0] = x_vec[x_vec.size()-1][0];
	    x0[1] = x_vec[x_vec.size()-1][1];
	    x0[2] = x_vec[x_vec.size()-1][2];

	    xmax = 0.0;
	    xmin = 0.0;

	    W  = 0.0;
	    i  = npoints-1;
	    rolling_t = 0.0;
	    while(rolling_t<=10*M_PI/wext)
	    {
	    	if(xmax<x_vec[i][1])
	    		xmax = x_vec[i][1];
	    	if(xmin>x_vec[i][1])
	    		xmin = x_vec[i][1];
	    	//W+=pow(V0/Rl - x_vec[i][2]/(CapEM(x_vec[i][1])*Rl),2)*Rl*(times[i]-times[i-1]);
	    	rolling_t+=times[i]-times[i-1];
	    	i--;
	    }

	    th_W.push_back(W);

	}
	return( dif(th_W, exp_W));
}

void export_curve( std::string export_filename, std::string waveform, double f_min, double f_max, double df )
{
	ofstream ofs(export_filename);
    ofstream ofs2(waveform);

   	// initial conditions
    state_type x0 = { 0.0 , 0.0, C0*V0};

    vector<state_type> x_vec;
	vector<double> times;

	double tmax;
	double dt;
	int npoints;

	double xmax;
	double xmin;

	double W;
	int i;
	double rolling_t;

   	for(double fex = f_min; fex<f_max; fex+=df)
    {
    	// external circular frequency
    	wext = 2*M_PI*fex;

    	// maximal time
    	tmax = 200.0*2*M_PI/wext;
    	dt = 1e-5*2*M_PI/wext;

	    // RK4 method
	    npoints = integrate_adaptive(controlled_stepper,oscillator,
	    	x0 , 0.0 , tmax , dt, push_back_state_and_time(x_vec, times));

	    x0[0] = x_vec[x_vec.size()-1][0];
	    x0[1] = x_vec[x_vec.size()-1][1];
	    x0[2] = x_vec[x_vec.size()-1][2];

	    xmax = 0.0;
	    xmin = 0.0;

	    W  = 0.0;
	    i  = npoints-1;
	    rolling_t = 0.0;
        ofs2<<"fw"<<fex<<std::endl;
	    while(rolling_t<=50*M_PI/wext)
	    {
            if(xmax<x_vec[i][1])
	    		xmax = x_vec[i][1];
	    	if(xmin>x_vec[i][1])
	    		xmin = x_vec[i][1];
            W+=pow(V0/Rl - x_vec[i][2]/(Cap(x_vec[i][1])*Rl),2)*Rl*(times[i]-times[i-1]);
	    	ofs2<<rolling_t<<","<<x_vec[i][1]<<","<<x_vec[i][0]<<","<<x_vec[i][2]<<std::endl;
            rolling_t+=times[i]-times[i-1];
	    	i--;
	    }

	    ofs<<fex<<","<<W/(fex*rolling_t)<<","<<0.5*(xmax-xmin)<<endl;

	}

	x0 = { 0.0 , 0.0, C0*V0};

   	for(double fex = f_max; fex>f_min; fex-=df)
    {
    	// external circular frequency
    	wext = 2*M_PI*fex;

    	// maximal time
    	tmax = 200.0*2*M_PI/wext;
    	dt = solv_delta*2*M_PI/wext;

	    // RK4 method
	    npoints = integrate_adaptive(controlled_stepper,oscillator,
	    	x0 , 0.0 , tmax , dt, push_back_state_and_time(x_vec, times));

	    x0[0] = x_vec[x_vec.size()-1][0];
	    x0[1] = x_vec[x_vec.size()-1][1];
	    x0[2] = x_vec[x_vec.size()-1][2];

	    xmax = 0.0;
	    xmin = 0.0;

	    W  = 0.0;
	    i  = npoints-1;
	    rolling_t = 0.0;
        ofs2<<"bw"<<fex<<std::endl;        
	    while(rolling_t<=50*M_PI/wext)
	    {
	    	if(xmax<x_vec[i][1])
	    		xmax = x_vec[i][1];
	    	if(xmin>x_vec[i][1])
	    		xmin = x_vec[i][1];
            //W+=x_vec[i][0]*Ftrans(x_vec[i][1],x_vec[i][2])*(times[i]-times[i-1]);
	    	W+=pow(V0/Rl - x_vec[i][2]/(Cap(x_vec[i][1])*Rl),2)*Rl*(times[i]-times[i-1]);
	    	ofs2<<rolling_t<<","<<x_vec[i][1]<<","<<x_vec[i][0]<<","<<x_vec[i][2]<<std::endl;
	    	rolling_t+=times[i]-times[i-1];
	    	i--;
	    }

	    ofs<<fex<<","<<W/(fex*rolling_t)<<","<<0.5*(xmax-xmin)<<endl;

	}
    ofs2.close();
	ofs.close();

}


int main(int argc, char **argv)
{
	Aext = 1.0*9.81;
	V0 = 21.0;
	export_curve("export_0V_1.0g.csv","wf_0V_1.0g.csv",20.0,350.0,2.0);

	V0 = 46.0;
	export_curve("export_25V_1.0g.csv","wf_25V_1.0g.csv",20.0,350.0,2.0);

	Aext = 0.5*9.81;
	V0 = 21.0;
	export_curve("export_0V_0.5g.csv","wf_0V_0.5g.csv",20.0,350.0,2.0);

	V0 = 46.0;
	export_curve("export_25V_0.5g.csv","wf_25V_0.5g.csv",20.0,350.0,2.0);

	Aext = 2.0*9.81;
	V0 = 21.0;
	export_curve("export_0V_2.0g.csv","wf_0V_2.0g.csv",20.0,450.0,2.0);

	V0 = 46.0;
	export_curve("export_25V_2.0g.csv","wf_25V_2.0g.csv",20.0,450.0,2.0);
}
