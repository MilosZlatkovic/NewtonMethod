#include <iostream>
#include <string>
#include <math.h>
#include <ginac/ginac.h> 
#include <vector>

using namespace std;
using namespace GiNaC;

static map<string, symbol> directory;

double norm(matrix a, int dim)
{
  double sum = 0.0;

  for(int i = 0; i < dim; i++)
    sum += ex_to<numeric>(a[i]).to_double() * ex_to<numeric>(a[i]).to_double();

  return sqrt(sum);
}


const symbol get_symbol(const string s)
{
  map<string, symbol>::iterator i = directory.find(s);
  if (i != directory.end())
    return i->second;
  else
    return directory.insert(make_pair(s, symbol(s))).first->second;
}


matrix calculateFunctionsMatrix(vector<ex> f, matrix x, int dim)
{
  int i;
  map<string, symbol>::iterator mapIter;
  vector<ex>::iterator vectorIter;
  
  exmap m;
  lst rez;

  for (i = 0, mapIter = directory.begin(); mapIter != directory.end() && i < dim; ++mapIter, i++) 
    m[get_symbol((mapIter -> first))] = ex_to<numeric>(x[i]).to_double();
 
  for(vectorIter = f.begin(); vectorIter != f.end(); ++vectorIter)
    rez.append((*vectorIter).subs(m));

  return matrix(dim, 1, rez);
}


matrix calculateJacobianMatrix(matrix a, matrix x, int dim)
{
  int i;
  map<string, symbol>::iterator mapIter;
  exmap m;
  lst rez;

  for (mapIter = directory.begin(), i = 0; mapIter != directory.end() && i< dim; ++mapIter, i++) 
    m[get_symbol((mapIter -> first))] = ex_to<numeric>(x[i]).to_double();
  
  for(i = 0; i < dim*dim; i++)
    rez.append(a[i].subs(m));	

  return matrix(dim, dim, rez);
}

int main()
{
  int max_iter, iter = 0;           
  double eps, temp;             
  vector<ex> f;        
  lst l, jacobian;

  //===================================================================================================== 
  // GiNaC provides no way to directly read an expression from a stream because you will usually want 
  // the user to be able to enter something like ‘2*x+sin(y)’ and have the ‘x’ and ‘y’ correspond to 
  // the symbols x and y you defined in your program and there is no way to specify the desired symbols 
  // to the >> stream input operator.
  
  // Set dimension;
   int dim  = 2;

  // Enter expressions;
  ex e1 = pow(get_symbol("x"),3) + pow(get_symbol("y"),3) - 1;
  ex e2 = pow(get_symbol("x"),2) + pow(get_symbol("y"),2) - 2*get_symbol("x");
  //=====================================================================================================
 
  f.push_back(e1);
  f.push_back(e2);

  cout << "Enter X0: " << endl;
 
  for(int i = 0; i < dim; i++)
    {
      cin >> temp;
      l.append(temp);
    }

  matrix X0(dim, 1, l);
  matrix X(dim, 1);
  matrix Xp(dim, 1);

  cout << "Enter the maximum number of iterations: " << endl;
  cin >> max_iter;

  cout << "Enter epsilon: " << endl;
  cin >> eps;
  
  for(int i = 0; i < dim; i++)
    {
      map<string, symbol>::iterator iter;
      for (iter = directory.begin(); iter != directory.end(); ++iter) 
	jacobian.append(f[i].diff(get_symbol(iter -> first)));	
    }

  matrix jacobianMatrix(dim, dim, jacobian);  
  
  X = X0;
  
  while(1)
    {
      iter++;
    
      if(iter >= max_iter)
	{
	  cout << "Method is divergent for given max number of iterations. " << endl;
	  exit(0);
	}

      Xp = X;

      X = Xp.sub( calculateJacobianMatrix(jacobianMatrix, Xp, dim).inverse().mul(calculateFunctionsMatrix(f,Xp, dim)));
      
      if(norm(X.sub(Xp), dim) < eps)
	break;   
    }

  cout << "Result: " << endl;
  cout << X << endl;
  
  return 0;
}
