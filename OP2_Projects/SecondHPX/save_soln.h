/*inline void save_soln(const double *q, double *qold){
  for (int n=0; n<4; n++) qold[n] = q[n];
}*/
#include <vector>

inline std::vector<double> save_soln(const double *q, double *qold){
	std::vector<double> result;

	for (int n=0; n<4; n++) qold[n] = q[n];

  	for (int i=0; i<4; ++i)
		result.push_back(qold[i]);

   	return result;
}


