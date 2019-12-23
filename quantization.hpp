#ifndef QUANTI
#define QUANTI




using namespace std;

//discretization from GATB     
// discretization scheme to store abundance values from 0 to 50000 on 8 bits
// with  5% error maximum
// from 0     to 70     :   step = 1              (70 buckets)
// from 70    to 100    :   step = 2              (15 buckets)
// from 100   to 500    :   step = 10             (40 buckets)
// from 500   to 1000   :   step = 20             (25 buckets)
// from 1000  to 5000   :   step = 100            (40 buckets)
// from 5000  to 10000  :   step = 200            (25 buckets)
// from 10000 to 50000  :   step = 1000           (40 buckets)
//
//to change discretization scheme, change the values in abundance_discretization below
MAX_ABUNDANCE_DISCRETE=65000; 

uint16_t abundance_at (uint8_t index, vector<uint16_t>& abudance_discretization)  
{
	if (index < abundance_discretization.size())
	{
		return floorf((abundance_discretization[index]  +  abundance_discretization[index+1])/2.0);
	}
	else 
	{
		return 0;
	}
}

vector<uint16_t> init_discretization_scheme()
{
	vector<uint8_t> abundance_discretization;
	abundance_discretization.resize(257);
	int total =0;
	abudance_discretization[0] = 0;
	int idx=1;
	for(int ii=1; ii<= 70; ii++,idx++ )
	{
		total += 1;
		abudance_discretization[idx] = total ;
	}
	
	for(int ii=1; ii<= 15; ii++,idx++ )
	{
		total += 2;
		abudance_discretization[idx] = total  ;
	}
	
	for(int ii=1; ii<= 40; ii++,idx++ )
	{
		total += 10;
		abudance_discretization[idx] = total  ;
	}
	for(int ii=1; ii<= 25; ii++,idx++ )
	{
		total += 20;
		abudance_discretization[idx] = total  ;
	}
	for(int ii=1; ii<= 40; ii++,idx++ )
	{
		total += 100;
		abudance_discretization[idx] = total  ;
	}
	for(int ii=1; ii<= 25; ii++,idx++ )
	{
		total += 200;
		abudance_discretization[idx] = total  ;
	}
	for(int ii=1; ii<= 40; ii++,idx++ )
	{
		total += 1000;
		abudance_discretization[idx] = total  ;
	}
	abudance_discretization[256] = total;
	return abundance_discretization;
}


int return_count_bin(uint16_t abundance, vector<uint16_t>& abudance_discretization)
{						
	int idx ;
	if (abundance >= MAX_ABUNDANCE_DISCRETE)
	{
		//~ _nb_abundances_above_precision++;
		//std::cout << "found abundance larger than discrete: " << abundance << std::endl;
		idx = abudance_discretization.size() -2 ;
	}
	else
	{
		//get first cell strictly greater than abundance
		vector<int>::iterator  up = upper_bound(abudance_discretization.begin(), abudance_discretization.end(), abundance);
		up--; // get previous cell
		idx = up- abudance_discretization.begin() ;
	}
	return idx;
}

#endif
