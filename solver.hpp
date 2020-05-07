#include <iostream>
#include <complex>

using namespace std;

namespace solver
{
    

    struct RealVariable
    {
        double n1,n2,n3;
        bool flag;
        RealVariable()
        {
            cout<<"new"<<endl;
            n1=0;
            n2=1;
            n3=0;
            flag=true;
        }

        

        RealVariable& operator * (RealVariable &other);
        friend RealVariable& operator * (RealVariable &other, double n);
        friend RealVariable& operator * (double n, RealVariable &other);
        

        RealVariable& operator - (RealVariable &other);
        friend RealVariable& operator - (RealVariable &other, double n);
        friend RealVariable& operator - (double n, RealVariable &other);
        
    
        RealVariable& operator ^ (double n);

        
        RealVariable& operator + (RealVariable &other);
        friend RealVariable& operator + (RealVariable &other, double n);
        friend RealVariable& operator + (double n, RealVariable &other);

        RealVariable& operator / (RealVariable &other);
        friend RealVariable& operator / (RealVariable &other, double n);
        friend RealVariable& operator / (double n, RealVariable &other);

        RealVariable& operator == (RealVariable &other);
        friend RealVariable& operator == (RealVariable &other, double n);
        friend RealVariable& operator == (double n, RealVariable &other);
      
    };
    double solve(RealVariable &other);
    double solve(double n);
   
    

    
    


    struct ComplexVariable
    {
        complex<double> n1, n2, n3;
        ComplexVariable()
        {
           n1 = complex<double>(0,0);
           n2 = complex<double>(1,0);
           n3 = complex<double>(0,0);
        }


        ComplexVariable& operator * (ComplexVariable &other);
        friend ComplexVariable& operator * (ComplexVariable &other, double n);
        friend ComplexVariable& operator * (double n, ComplexVariable &other);
        friend ComplexVariable& operator * (complex<double> c, ComplexVariable &other);
        friend ComplexVariable& operator * (ComplexVariable &other ,complex<double> c);

        ComplexVariable& operator - (ComplexVariable &other);
        friend ComplexVariable& operator - (ComplexVariable &other, double n);
        friend ComplexVariable& operator - (double n, ComplexVariable &other);
        friend ComplexVariable& operator - (complex<double> c, ComplexVariable &other);
        friend ComplexVariable& operator - (ComplexVariable &other ,complex<double> c);

        ComplexVariable& operator == (ComplexVariable &other);
        friend ComplexVariable& operator == (ComplexVariable &other, double n);
        friend ComplexVariable& operator == (double n, ComplexVariable &other);
        friend ComplexVariable& operator == (complex<double> c, ComplexVariable &other);
        friend ComplexVariable& operator == (ComplexVariable &other ,complex<double> c);

        ComplexVariable& operator ^ (double n);

        ComplexVariable& operator + (ComplexVariable &other);
        friend ComplexVariable& operator + (ComplexVariable &other, double n);
        friend ComplexVariable& operator + (double n, ComplexVariable &other);
        friend ComplexVariable& operator + (complex<double> c, ComplexVariable &other);
        friend ComplexVariable& operator + (ComplexVariable &other ,complex<double> c);

        ComplexVariable& operator / (ComplexVariable &other);
        friend ComplexVariable& operator / (ComplexVariable &other, double n);
        friend ComplexVariable& operator / (double n, ComplexVariable &other);
        friend ComplexVariable& operator / (complex<double> c, ComplexVariable &other);
        friend ComplexVariable& operator / (ComplexVariable &other ,complex<double> c);

    };
    complex<double> solve(ComplexVariable &c);
    

    

    ComplexVariable& operator - (double n, ComplexVariable &other);

    ComplexVariable& operator + (double n ,ComplexVariable &other);
    ComplexVariable& operator + (ComplexVariable& other ,complex<double> other2);
    ComplexVariable& operator + (complex<double> &other);
    
    
   

};