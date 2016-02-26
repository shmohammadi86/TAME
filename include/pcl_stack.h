
#ifndef _STACK_PCL_
#define _STACK_PCL_

//#include "utils.h"


// MUST call the init function to allocate memory for the queue
// size has be a power of 2 always
 
template <class data_type>
class stack_pcl
{
public:
	stack_pcl():m_data(NULL){ }
	~stack_pcl()
	{
		if(m_data != NULL)
		{	
			free(m_data);
			m_data = NULL;
		}
	}

	void init(unsigned int size)
	{
		m_data = (data_type*) malloc(sizeof(data_type) * size);
		m_top = size - 1;
		m_count = 0;
		
		m_stack_size = size;		
		m_stack_size_1 = size-1;		
	}

	void reset()
	{
		m_top = m_stack_size - 1;
		m_count = 0;
	}	

	void push(data_type& x)
	{
		//if (m_count >= m_stack_size)
			//cout << "Warning: stack overflow push x " << endl;
	        //else 
		{
			m_top = (m_top + 1) & (m_stack_size_1);
			m_data[m_top] = x;			
			m_count = m_count + 1;
		}	
	}

	data_type pop()
	{
		data_type value;

		//if (m_count <= 0) 
			//cout << "Warning: empty stack pop." << endl;
		//else 
		{
			value = m_data[m_top--];
			m_top = (m_top & m_stack_size_1);
			m_count = m_count - 1;			
		}
		
		return value;
	}

	bool empty()
	{
		return (m_count == 0);	
	}

	int size()
	{
		return m_count;	
	}

public:
	data_type*	m_data;		
	int 		m_top;                    
	int 		m_count;                    

	unsigned int    m_stack_size; // this has to be 2^x
	unsigned int    m_stack_size_1; // this has to be 2^x-1 
};

#endif
