#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#define UPPER_LIMIT 1000000000
using namespace std;

int main(int argc, char** argv)
{
	string infilename, outfilename;
	
	if (argc > 2)
	{
        	infilename = argv[1];
		outfilename = argv[2];
	}
    	else
    	{
        	cout << "./mtx2bin input_filename output_filename!!!\n";
        	return -1;
    	}

	ifstream in(infilename.c_str());
        if(!in)
        {
                cout << infilename <<" not Found!" << endl;
                return 0;
        }

        istringstream in2;
        string line="";

        getline(in,line);
        while(line.size() > 0 && line[0] == '%')
                getline(in,line);

        int m, n, v_count;
        int64_t e_count;

        in2.str(line);
        in2 >> m >> n >> e_count;

        if(m !=n)
        {
                cout << " m != n" << endl;
                return false;
        }

        v_count = m;
        vector <int> R;
	vector <double> W;	
	R.reserve(e_count + e_count);
	W.reserve(e_count);	

        int64_t entry_counter = 0;
        int target;
        int source;
        int64_t exact_e_count = 0;//, 
	double	weight;

        while(!in.eof() && entry_counter < e_count)
        {
                getline(in, line);
                entry_counter++;

                if(line!="")
                {
                        in2.clear();
                        in2.str(line);

                        in2 >> source >> target >> weight;

                        if (source <= target)
                        {
                                //cout << target << " " << source << endl;
                                //cout << "found a nnz in the upper triangle (self loop) and ignored" << endl;
                                continue;
                        }

			R.push_back(source);
			R.push_back(target);
			W.push_back(weight);	
                }
        }

        in.close();


	cout << "writing" << endl;
	ofstream filew (outfilename.c_str(), ios::out|ios::binary);
  	if (filew.is_open())
  	{
		filew.write((char*)&v_count, sizeof(int));
		e_count = R.size()/2;
		filew.write((char*)&e_count, sizeof(int64_t));

		cout << "writing " << v_count << " " << e_count << endl;	

		for(int64_t i = 0; i < 2; i++)
		{	
			for(int64_t j = 0; j < 2; j++)
				cout << R[i * 2 + j] << " ";
			cout << W[i];
			cout << endl;	
		}
 		
		for(int64_t i = e_count - 2; i < e_count; i++)
		{	
			for(int64_t j = 0; j < 2; j++)
				cout << R[i * (int64_t)2 + j] << " ";
			cout << W[i];	
			cout << endl;	
		}
		int64_t size = R.size();    	   	
		int64_t w_size = W.size();

		filew.write((char*)&R[0], sizeof(int) * size);
		filew.write((char*)&W[0], sizeof(double) * w_size);
    		filew.close();
  	}
  	else 
		cout << "Unable to open file";

	cout << "DONE" << endl;
	R.clear();

	return 0;
}
