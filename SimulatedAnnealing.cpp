// Simulated Annealing
// Joshua D. John
// November 2021
// 
// Compile it with:
//   g++ -o SimulatedAnnealing SimulatedAnnealing.cpp -lboost_iostreams -lboost_system -lboost_filesystem

#include <iostream>
#include <fstream>
//#include <string>
#include <iomanip>
#include <map>
#include <utility>
#include <vector>
#include <cmath>

#include "gnuplot-iostream.h"

class Annealer
{
	private:
		typedef std::pair<float, float> pair;
		float surface_size{2};
		float grid_size{0.01};
		size_t n_rows = (int)surface_size/grid_size;
		size_t n_cols= (int)surface_size/grid_size;
		size_t n_size = n_rows * n_cols;
			
		int n = 10; // Number of cycles
		int m = 100; // Number of trials per cycle
		int na = 0; // Number of accepted solutions
		float p1 = 0.7; // Probability of accepting worse solution at the start
		float p50 = 0.001; // Probability of accepting worse solution at the end
		float t1 = -1.0/log(p1); // Initial temperature
		float t50 = -1.0/log(p50); // Final temperature
		float num = (t50/t1);
		float denom = (1.0/(n-1.0));
		float frac = pow(num,denom); // Fractional reduction every cycle
		
		/* float x_start{0.8};
		float y_start{-0.5};
		
		float* surface;
		float** surface_grid;
		float* x_grid;
		float* y_grid;
		*/
		
		std::vector<float> surface;
		std::vector<float> x_grid;
		std::vector<float> y_grid;
		std::vector<pair> x_list;
		std::vector<float> fs;
		std::vector<std::tuple<double, double, double>> optim_path;
		
		
		
		// std::vector<pair> x_start = {std::make_pair(0.8, -0.5)};
 		
		float snap(float value);
		double round_up(double value, int decimal_places);
	
	public:
		Annealer();
		~Annealer();
		
		void surface_one();
		void plot();
		void optimize();
};

Annealer::Annealer()
{
	/*this->surface = new float[this->n_size];
	this->x_grid = new float[this->n_rows];
	this->y_grid = new float[this->n_cols];
	
	int i = 0;
	for (float x = -0.5*this->surface_size; x < 0.5*this->surface_size; x+=this->grid_size)
	{
		this->x_grid[i] = x;
		this->y_grid[i] = x;
		i++;
	}*/
	
	for (float x = -0.5*this->surface_size; x < 0.5*this->surface_size; x+=this->grid_size)
	{
		this->x_grid.push_back(x);
		this->y_grid.push_back(x);
	}	
}

Annealer::~Annealer()
{

}

float Annealer::snap(float value)
{
    // Added std::abs to give correct behaviour for negative values
    return value - std::abs(std::fmod(value, this->grid_size));
}

double Annealer::round_up(double value, int decimal_places) 
{
    const double multiplier = std::pow(10.0, decimal_places);
    return std::ceil(value * multiplier) / multiplier;
}

void Annealer::surface_one()
{	
	for (size_t x {0}; x < this->n_rows; x++)
	{
		for (size_t y {0}; y < this->n_cols; y++)
		{
			int xy = (x*this->n_rows) + y;

			float xx = this->x_grid[x];
			float yy = this->y_grid[y];
             		float z = 0.2 + xx*xx + yy*yy - 0.1*cos(6.0*M_PI*xx) - 0.1*cos(6.0*M_PI*yy);
			this->surface.push_back(z);
		}
	}
}

void Annealer::optimize()
{
	std::map<pair, float> surface_map;
	auto x_start = std::make_pair(0.8, -0.5);
	
	for (size_t x {0}; x < this->n_rows; x++)
	{
		for (size_t y {0}; y < this->n_cols; y++)
		{
			int xy = (x*n_rows) + y;
			
			float xx = this->round_up(this->x_grid[x],2);
			float yy = this->round_up(this->y_grid[y],2);
			float zz = this->surface[xy];
			surface_map.insert({ std::make_pair(xx, yy), zz });
		}
	}
	
	/*for (const auto &entry: surface_map)
	{
		auto key_pair = entry.first;
		std::cout << "{" << key_pair.first << "," << key_pair.second << "}, "
		          << entry.second << std::endl;
	}
	
	auto search = surface_map.find(x_start);
    	if (search != surface_map.end()) 
    	{
    		auto key_pair = search->first;
        	std::cout << "Found " << "{" << key_pair.first << "," << key_pair.second << "}, "
        		 << search->second << '\n';
    	} 
    	else 
    	{
        	std::cout << "Not found\n";
    	}
    	
    	std::cout << '(' << std::get<0>(x_start) << ", " << std::get<1>(x_start) << ")\n";
    	x_start = std::make_pair(0.4, 0.1);
    	std::cout << '(' << std::get<0>(x_start) << ", " << std::get<1>(x_start) << ")\n";
    	*/
    	
    	bool accept;
    	
    	this->x_list.push_back(x_start);
    	auto x_init = x_start;
    	
    	this->na += 1.0;
    	auto x_current = x_list[0];
    	auto f = surface_map.find(x_init);
    	float f_current = f->second;
    	
    	this->fs.push_back(f_current);
    	
    	float t = this->t1;
    	float DeltaE_avg = 0.0;
    	float p;
    	
    	for (int i{0}; i < n; i++)
    	{
    		// std::cout << "Cycle: " << i << " - Temperature: " << t << std::endl;
    		for (int j{0}; j < m; j++)
    		{
    			double xs1 = std::get<0>(x_init) +  (double)rand()/RAND_MAX - 0.5;
    			double xs2 = std::get<1>(x_init) +  (double)rand()/RAND_MAX - 0.5;
    			//std::cout << "xs-1: (" << xs1 << "," << xs2 << ")\n";
    			//clip between -1,1
    			xs1 = std::max(std::min(xs1, 0.99), -0.99);
        		xs2 = std::max(std::min(xs2, 0.99), -0.99);
        		//std::cout << "xs-clip: (" << xs1 << "," << xs2 << ")\n";
        		xs1 = this->round_up(this->snap(xs1), 2);
        		xs2 = this->round_up(this->snap(xs2), 2);
        		//std::cout << "xs-snap: (" << xs1 << "," << xs2 << ")\n";
        		
        		x_init = std::make_pair(xs1, xs2);
        		auto ff = surface_map.find(x_init);
        		float f_x = ff->second; 
        		float DeltaE = fabs(f_x - f_current);
        		//std::cout << "DeltaE: " << DeltaE << std::endl; 
        		if (f_x > f_current)
        		{
        			// Initialize DeltaE_avg if a worse solution was found on the first iteration
        			if (i==1 && j==1)
        			{
                			DeltaE_avg = DeltaE;
                		}
                		
                		// objective function is worse, generate probability of acceptance
                		p = exp(-DeltaE/(DeltaE_avg * t));
                		
        		}
        		else
        		{
        			if (((float)rand()/RAND_MAX) < p)
        			{
                			// accept the worse solution
                			accept = true;
                		}
            			else
                		{
                			// don't accept the worse solution
                			accept = false;
                		}	
                		// objective function is lower, automatically accept
            			accept = true;
        		}
        		
        		if (accept == true)
        		{
        			x_current = x_init;
        			auto f2 = surface_map.find(x_current);
        			float f_xc = f2->second; 
        			f_current = f_xc;
        			
        			// increment number of accepted solutions
            			na += 1.0;
            			// update DeltaE_avg
            			DeltaE_avg = (DeltaE_avg * (na - 1.0) +  DeltaE) / na;
        		}
    			   		
    		}
    		
    		// Record the best x values at the end of every cycle
    		this->x_list.push_back(x_current);
    		this->fs.push_back(f_current);
   		this->optim_path.push_back(std::make_tuple(std::get<0>(x_current), std::get<1>(x_current), f_current));
   		// Lower the temperature for next cycle
    		t = frac * t;
    	}
    	std::cout << "Best solution: " << "(" << std::get<0>(x_current) << ", " << std::get<1>(x_current) << ")\n";
	std::cout << "Best objective: " << f_current << std::endl;	
}

void Annealer::plot()
{
	Gnuplot gp;
	std::vector<std::tuple<double, double, double>> surfaceGridPnts;
	std::vector<std::tuple<double, double, double, double, double, double>> optim_vector;
	
	for (int i {0}; i < (this->optim_path.size() - 1); i++)
	{
		double x1 = std::get<0>(this->optim_path[i+1]);
		double y1 = std::get<1>(this->optim_path[i+1]);
		double z1 = std::get<2>(this->optim_path[i+1]);
		double x2 = std::get<0>(this->optim_path[i]);
		double y2 = std::get<1>(this->optim_path[i]);
		double z2 = std::get<2>(this->optim_path[i]);
		
		optim_vector.push_back(std::make_tuple(x1, y1, z1, x2, y2, z2));
	}
	for (size_t x{0}; x < this->n_rows; x++)
	{
		for (size_t y{0}; y < this->n_cols; y++)
		{
			int xy = (x*n_rows) + y;
			
			float xx = this->x_grid[x];
			float yy = this->y_grid[y];
			float zz = this->surface[xy];
			surfaceGridPnts.push_back(std::make_tuple(xx, yy, zz));
		}
	}
	
	/*std::vector<std::tuple<double, double, double>> optim_path;
	for (size_t i{0}; i < n; i++)
	{
		float xp = std::get<0>(this->x_list[i]);
		float yp = std::get<0>(this->x_list[i]);
		float zp =  this->fs[i];
		optim_path.push_back(xp, yp, zp)
	}*/

	/*for(const auto& value: this->optim_path) 
	{
    		std::cout << "(" << std::get<0>(value) << "," << std::get<1>(value)  << ")" << "\n";
	}*/
	
	/*gp << "set samples 30\n";
	// gp << "set isosamples 30 \n";
	// gp << "set pm3d interpolate 1,1\n";
	// gp << "splot [-2:2][-2:2] " << "++"  << "u 1:2:(f($1, $2)) w l\n";
	gp << "set view map\n";
	gp << "set size ratio -1\n";
	gp << "set surface\n";
	gp << "set dgrid3d 10,10 qnorm 2\n";
	gp << "set contour base\n";
	gp << "set cntrlabel format '%4.3g' font ',7' start 5 interval 10\n";
	gp << "set cntrparam bspline\n";
	gp << "set key outside\n";
	//gp << "set style data lines\n";
	gp << "set isosamples 30\n";
	gp << "unset colorbox\n";
	// gp << "set cbrange [0:7000]\n";
	gp << "set cntrparam levels 10\n";
	// gp << "set palette model RGB defined ( 0 'white', 1 'black' )\n";
	// gp << "set style line 1 lc rgb '#4169E1' pt 7 ps 2\n";
	gp << "splot" << gp.file1d(surfaceGridPnts) << "with lines title 'Optimize', " 
		<< gp.file1d(this->optim_path) << "using 1:2:(0) with points linecolor rgb ' white'notitle\n";*/
	
	gp << "set xrange [-1:1]\nset yrange [-1:1]\nset zrange [0:2]\n";
	//gp << "set hidden3d\n";
	gp << "set pm3d interpolate 1,1\n";
	/*gp << "splot" << gp.file1d(surfaceGridPnts) << "with lines lc rgb 'blue' title 'Simulated Annealing'," << gp.file1d(this->optim_path) << "using 1:2:3 with linespoints linecolor rgb 'red' notitle\n";*/
	/*gp << "splot" << gp.file1d(surfaceGridPnts) << "with lines lc rgb 'blue' title 'Simulated Annealing'," << gp.file1d(optim_vector) << "using 1:2:3:(1):(1):(1) with vector linecolor rgb 'red' notitle\n";*/
	gp << "splot" << gp.file1d(surfaceGridPnts) << "with lines lc rgb 'pink' title 'Simulated Annealing'," << gp.file1d(this->optim_path) << "pt 7 title 'Iteration number',"  << gp.file1d(this->optim_path) << "using 1:2:3:($0+1) with labels offset 1 notitle\n";
	
}

int main()
{
	Annealer al{};
	al.surface_one();
	al.optimize();
	al.plot();
	
	return 0;
}
