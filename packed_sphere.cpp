#include <iostream>
#include <math.h>

using namespace std;

int main()
{
	int side = 4;
	for(int i = 0; i < side; ++i)
	{
		for(int j = 0; j < side; ++j)
		{
			for(int k = 0; k < side; ++k)
			{
				double x = 2 * i + j%2 + k%2;
				double y = sqrt(3.0) * (j + ((1.0/3.0) * (k%2)));
				double z = (2 * sqrt(6))/3.0 * k;
				cout<<x<<" "<<y<<" "<<z<<endl;

			}
		}
	}

}
