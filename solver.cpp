#include "solver.hpp"
#include <iostream>
#include <stdexcept>
#include <complex>
#include <math.h>

using  namespace std;
using namespace solver;


  //RealVariable tmp; 
  double num=0, cx=0;
    
    double solver::solve(RealVariable &other)
    {
        cout<<"in"<<endl;
        return num;
    }
    double solver::solve(double n)
    {
        
        return 0;
    }



    RealVariable& solver::operator * (RealVariable &other, double n)
    {   
        cout<<"*"<<endl;
        RealVariable *tmp = new RealVariable;
        //n = num*n;

        tmp->n1= other.n1 ;
        tmp->n2 = other.n2*n;
        tmp->n3= other.n3 ;
        cout<<"1)*)"<<tmp->n1<<" , "<<tmp->n2<<" , "<<tmp->n3<<endl;
        RealVariable &tmp2=*tmp;
        return tmp2;
        
    }   
    RealVariable& RealVariable::operator * (RealVariable &other)
    {   
        cout<<"*"<<endl;
        try{
            if(this->flag == false)
               throw std::invalid_argument("exception");

        }catch(const std::exception& e)
		{
			std::cerr << e.what() << '\n';
			throw;
		}
        
        RealVariable *tmp = new RealVariable;
        tmp->n1=other.n1+1;
        tmp->n2=0;
        tmp->n3=other.n3;            
        tmp->flag=false;

        RealVariable &tmp2=*tmp;
        cout<<"2)*)"<<tmp->n1<<" , "<<tmp->n2<<" , "<<tmp->n3<<endl;
        return tmp2;
        
        
    }
    RealVariable& solver::operator * (double n, RealVariable &other)
    {  
        cout<<"*"<<endl;
        cout<<other.n2*n<<endl;
        RealVariable *tmp = new RealVariable;
        tmp->n1 = other.n1*n;
        tmp->n2 = other.n2*n;
        tmp->n3 = other.n3*n;
        num=n;
        cout<<"3)*)"<<tmp->n1<<" , "<<tmp->n2<<" , "<<tmp->n3<<endl;
        RealVariable &tmp2=*tmp;
        return tmp2;
    }


    RealVariable& solver::operator - (RealVariable &other, double n)
    {   
        cout<<"-"<<endl;
        RealVariable *tmp = new RealVariable;
        
        tmp->n1=other.n1;
        tmp->n2=other.n2;
        tmp->n3=other.n3-n;

        cout<<"1)-)"<<tmp->n1<<" , "<<tmp->n2<<" , "<<tmp->n3<<endl;
    
        RealVariable &tmp2=*tmp;
        return tmp2;
        
    }   
    RealVariable& RealVariable::operator - (RealVariable &other)
    {   
        cout<<"-"<<endl;
        RealVariable *tmp = new RealVariable;
        tmp->n1=this->n1-other.n1;
        tmp->n2=this->n2-other.n2;
        tmp->n3=this->n3-other.n3;
        cout<<"2)-)"<<tmp->n1<<" , "<<tmp->n2<<" , "<<tmp->n3<<endl;
        RealVariable &tmp2=*tmp;
        return tmp2;
        
    }
    RealVariable& solver::operator - (double n, RealVariable &other)
    {   
        cout<<"-"<<endl;
        RealVariable *tmp = new RealVariable;
        
        tmp->n1 = -other.n1;
        tmp->n2 = -other.n2;
        tmp->n3 = n-other.n3;
        cout<<"3)-)"<<tmp->n1<<" , "<<tmp->n2<<" , "<<tmp->n3<<endl;
        num=1;
        RealVariable &tmp2=*tmp;
        return tmp2;
    }
    

   
    RealVariable& RealVariable::operator ^ (double n)
    {   
        try{
            if(n != 2)
                throw std::invalid_argument("exception");

        RealVariable *tmp = new RealVariable;
        tmp->n1=this->n1+1;
        tmp->n2=0;
        tmp->n3=this->n3;            

        RealVariable &tmp2=*tmp;
        cout<<"2)*)"<<tmp->n1<<" , "<<tmp->n2<<" , "<<tmp->n3<<endl;
        return tmp2;

        }catch(const std::exception& e)
		{
			std::cerr << e.what() << '\n';
			throw;
		}
       
    }



    RealVariable& solver::operator + (RealVariable &other, double n)
    {   
        cout<<"+"<<endl;
        RealVariable *tmp = new RealVariable;
        
        tmp->n1=other.n1;
        tmp->n2=other.n2;
        tmp->n3=other.n3+n;

        cout<<"1)+)"<<tmp->n1<<" , "<<tmp->n2<<" , "<<tmp->n3<<endl;
    
        RealVariable &tmp2=*tmp;
        return tmp2;
        
        
    }   
    RealVariable& RealVariable::operator + (RealVariable &other)
    {   
        cout<<"+"<<endl;
        RealVariable *tmp = new RealVariable;
        tmp->n1=this->n1+other.n1;
        tmp->n2=this->n2+other.n2;
        tmp->n3=this->n3+other.n3;
        cout<<"2)+)"<<tmp->n1<<" , "<<tmp->n2<<" , "<<tmp->n3<<endl;
        RealVariable &tmp2=*tmp;
        return tmp2;
        
    }
    RealVariable& solver::operator + (double n, RealVariable &other)
    {   
        cout<<"+"<<endl;
        RealVariable *tmp = new RealVariable;
        
        tmp->n1 = other.n1;
        tmp->n2 = other.n2;
        tmp->n3 = n+other.n3;
        cout<<"3)-)"<<tmp->n1<<" , "<<tmp->n2<<" , "<<tmp->n3<<endl;
    
        RealVariable &tmp2=*tmp;
        return tmp2;
    }



    RealVariable& solver::operator / (RealVariable &other, double n)
    {   
        cout<<"/"<<endl;
        RealVariable *tmp = new RealVariable;
        //n = num*n;
        if(n == 0 )
            throw std::invalid_argument("exception");
        tmp->n1 = other.n1 / n ;
        tmp->n2 = other.n2 / n;
        tmp->n3 = other.n3 / n;
        cout<<"1)/)"<<tmp->n1<<" , "<<tmp->n2<<" , "<<tmp->n3<<endl;
        RealVariable &tmp2=*tmp;
        return tmp2;
    
    }   
    RealVariable& RealVariable::operator / (RealVariable &other)
    {   
        cout<<"/"<<endl;
        RealVariable *tmp = new RealVariable;
        if(other.n1 != 0 && other.n2 != 0 && other.n3 !=0)
        {
            
            tmp->n1=this->n1/other.n1;
            tmp->n2=this->n2/other.n2;
            tmp->n3=this->n3/other.n3;

        }   
        else if(other.n1 == 0 && other.n2 == 0 && other.n3 != 0)
        {
            
            tmp->n1=this->n1;
            tmp->n2=this->n2;
            tmp->n3=this->n3/other.n3;

        }
        else if(other.n1 == 0 && other.n2 != 0 && other.n3 == 0)
        {
            
            tmp->n1=this->n1;
            tmp->n2=this->n2/other.n2;
            tmp->n3=this->n3;

        }
        else if(other.n1 == 0 && other.n2 != 0 && other.n3 != 0)
        {
            
            tmp->n1=this->n1;
            tmp->n2=this->n2/other.n2;
            tmp->n3=this->n3/other.n3;
        }
        else if(other.n1 != 0 && other.n2 == 0 && other.n3 == 0)
        {
           
            tmp->n1=this->n1/other.n1;
            tmp->n2=this->n2;
            tmp->n3=this->n3;
        }
        else if(other.n1 != 0 && other.n2 == 0 && other.n3 != 0)
        {
            
            tmp->n1=this->n1/other.n1;
            tmp->n2=this->n2;
            tmp->n3=this->n3/other.n3;
        }
        else if(other.n1 != 0 && other.n2 != 0 && other.n3 == 0)
        {
            
            tmp->n1=this->n1/other.n1;
            tmp->n2=this->n2/other.n2;
            tmp->n3=this->n3;
        }
        else
        {
            throw std::invalid_argument("exception");
        }
        
        cout<<"2)/)"<<tmp->n1<<" , "<<tmp->n2<<" , "<<tmp->n3<<endl;
       
        

        RealVariable &tmp2=*tmp;
        return tmp2;
        
    }
    RealVariable& solver::operator / (double n, RealVariable &other)
    {   
        cout<<"/"<<endl;
        cout<<other.n2*n<<endl;
        RealVariable *tmp = new RealVariable;
        tmp->n1 = n / other.n1;
        tmp->n2 = n / other.n2;
        tmp->n3 = n / other.n3;
        
        cout<<"3)/)"<<tmp->n1<<" , "<<tmp->n2<<" , "<<tmp->n3<<endl;
        RealVariable &tmp2=*tmp;
        return tmp2;
    }



    RealVariable& solver::operator == (RealVariable &other, double n)
    {   
        cout <<"=="<<endl;
        double a=0,b=0,c=0;
        RealVariable *tmp = new RealVariable;
        tmp->n1=other.n1;
        tmp->n2=other.n2;
        tmp->n3=other.n3-n;

        if(tmp->n1 != 0 )
         {
             num=0;
             a=tmp->n1;
             b=tmp->n2;
             c=tmp->n3;

            num=(-b+sqrt(b*b-4*a*c))/(2*a);
           
         }  
         else
         {
             num =0;
             b=tmp->n2;
             c=tmp->n3;

            if(b==0)
                throw std::invalid_argument("exception");   
            else if(c==0)
                num = 0;
            else   
                num = -(c/b);
             
             
         }
          
        return other; 
    }
    RealVariable& RealVariable::operator == (RealVariable &other)
    {
        cout <<"=="<<endl;
        double a=0,b=0,c=0;
        RealVariable *tmp = new RealVariable;
        tmp->n1=this->n1-other.n1;
        tmp->n2=this->n2-other.n2;
        tmp->n3=this->n3-other.n3;
        
        if(tmp->n1 != 0 )
         {
            
             num=0;
             a=tmp->n1;
             b=tmp->n2;
             c=tmp->n3;

            num=(-b+sqrt(b*b-4*a*c))/(2*a);
           
         }  
         else
         {
             
             num =0;
             b=tmp->n2;
             c=tmp->n3;

            if(b==0)
                throw std::invalid_argument("exception");   
            else if(c==0)
                num = 0;
            else   
                num = -(c/b);
             
         }
          
        
        return *this;
    }
    RealVariable& solver::operator == (double n, RealVariable &other)
    {
        cout <<"=="<<endl;
        other == n;
        
        return other; 
    }

    

    complex<double> solver::solve(ComplexVariable &other)
    {
        cout<<other.n1<<" , "<<other.n2<<" , "<<other.n3<<endl;
        complex<double> a,b,c,num1,num2;
        a=other.n1;
        b=other.n2;
        c=other.n3;

        if(a == complex(0.,0.) && b == complex(0.,0.))
        {
            throw std::invalid_argument("exception");   
        }
        else if(a == complex(0.,0.) && b != complex(0.,0.))
        {
            
            num1=-(c/b);
          
            return num1;
        }
        else
        { 
            num2= (other.n2*other.n2) - (4.*other.n1*other.n3);

            if(num2.real() != 0.0 && num2.imag() == 0)
            {
                num1 = (-other.n2 + sqrt(num2)) / (2.*other.n1);
                return num1;
            }

            if(num2 == 0.0)
            {
                num1 = -(other.n2/(2.*other.n1));
                return num1;
            }
        }
          return num1;     
    }         
    


    ComplexVariable& ComplexVariable::operator * (ComplexVariable &other)
    {
        ComplexVariable *tmp = new ComplexVariable;
  
    
        tmp->n1=(this->n1*other.n3)+(other.n1*this->n3)+(this->n2*other.n2);
        tmp->n2=(this->n2*other.n3)+(other.n2*this->n3);
        tmp->n3=this->n3*other.n3;

        ComplexVariable &tmp2=*tmp;

        return tmp2;
    }
    ComplexVariable& solver::operator * (ComplexVariable &other, double n)
    {
        ComplexVariable *tmp = new ComplexVariable;

        tmp->n1 = other.n1 * n;
        tmp->n2 = other.n2 * n;
        tmp->n3 = other.n3 * n;

        ComplexVariable &tmp2=*tmp; 
        return tmp2;
    }
    ComplexVariable& solver::operator * (double n, ComplexVariable &other)
    {
        ComplexVariable *tmp = new ComplexVariable;

        tmp->n1 = n * other.n1;
        tmp->n2 = n * other.n2;
        tmp->n3 = n * other.n3;
        cout<<tmp->n1<<" , "<<tmp->n2<<" , "<<tmp->n3<<endl;
        ComplexVariable &tmp2=*tmp;

        return tmp2;
    }
    ComplexVariable& operator*(complex<double> c, ComplexVariable &other)
    {
        ComplexVariable *tmp = new ComplexVariable;
    
        tmp->n1=other.n1*num;
        tmp->n2=other.n2*num;
        tmp->n3=other.n3*num;

        ComplexVariable &tmp2=*tmp;
    
        return tmp2;
    }
    ComplexVariable& operator*(ComplexVariable &other ,complex<double> c)
    {
        ComplexVariable *tmp = new ComplexVariable;
    
        tmp->n1= num * other.n1;
        tmp->n2= num * other.n2;
        tmp->n3= num * other.n3;

        ComplexVariable &tmp2=*tmp;

        return tmp2;
    }



    ComplexVariable& ComplexVariable::operator - (ComplexVariable &other)
    {
        ComplexVariable *tmp = new ComplexVariable;

        tmp->n1=this->n1-other.n1;
        tmp->n2=this->n2-other.n2;
        tmp->n3=this->n3-other.n3;

        ComplexVariable &tmp2=*tmp;
        
        return tmp2;
    }
    ComplexVariable& solver::operator - (ComplexVariable &other, double n)
    {
        ComplexVariable *tmp = new ComplexVariable;

        tmp->n1=other.n1;
        tmp->n2=other.n2;
        tmp->n3=other.n3-n;
        cout<<tmp->n1<<" , "<<tmp->n2<<" , "<<tmp->n3<<endl;
        ComplexVariable &tmp2=*tmp;

        return tmp2;
    }
    ComplexVariable& solver::operator - (double n, ComplexVariable &other)
    {
        ComplexVariable *tmp = new ComplexVariable;
        
        tmp->n1= -other.n1;
        tmp->n2= -other.n2;
        tmp->n3= n-other.n3;

        ComplexVariable& tmp2=*tmp;

        return tmp2;
    }
    ComplexVariable& solver::operator - (complex<double> c, ComplexVariable &other)
    {
        ComplexVariable *tmp = new ComplexVariable;

        tmp->n1= -other.n1;
        tmp->n2= -other.n2;
        tmp->n3= c-other.n3;

        ComplexVariable &tmp2=*tmp;

        return tmp2;
    }
    ComplexVariable& solver::operator - (ComplexVariable &other, complex<double> c)
    {
        ComplexVariable *tmp = new ComplexVariable;

        tmp->n1=other.n1;
        tmp->n2=other.n2;
        tmp->n3=other.n3-c;

        ComplexVariable &tmp2= *tmp;

        return tmp2;

    }



    ComplexVariable& ComplexVariable::operator == (ComplexVariable &other)
    {
        return *this-other;
    }
    ComplexVariable& solver::operator == (ComplexVariable &other, double n)
    {
        return other-n;
    }
    ComplexVariable& solver::operator == (double n, ComplexVariable &other)
    {
        return n-other;
    }
    ComplexVariable& solver::operator == (ComplexVariable &other, complex<double> c)
    {
        return other-c;
    }
    ComplexVariable& solver::operator == (complex<double> c, ComplexVariable &other)
    {
        return c-other;
    }
    



    ComplexVariable& ComplexVariable::operator ^ (double n)
    {
        try{
            if(n != 2)
                throw std::invalid_argument("exception");

            ComplexVariable *tmp = new ComplexVariable;

            tmp->n1 = 1;
            tmp->n2 = 0;
            tmp->n3 = 0;

            ComplexVariable &tmp2 =*tmp;

            return tmp2;

        }catch(const std::exception& e)
		{
			std::cerr << e.what() << '\n';
			throw;
		}
       
        
    }



    ComplexVariable& ComplexVariable::operator + (ComplexVariable &other)
    {
        ComplexVariable *tmp = new ComplexVariable;

        tmp->n1=this->n1 + other.n1;
        tmp->n2=this->n2 + other.n2;
        tmp->n3=this->n3 + other.n3;

        ComplexVariable &tmp2=*tmp;

        return tmp2;
    }
    ComplexVariable& solver::operator + (ComplexVariable &other, double n)
    {
        ComplexVariable *tmp = new ComplexVariable;
        tmp->n1=other.n1;
        tmp->n2=other.n2;
        tmp->n3=other.n3 + n;

        ComplexVariable &tmp2=*tmp;

        return tmp2;
    }
    ComplexVariable& solver::operator + (double n ,ComplexVariable &other)
    {
        ComplexVariable *tmp = new ComplexVariable;

        tmp->n1=other.n1;
        tmp->n2=other.n2;
        tmp->n3=n + other.n3;

        ComplexVariable &tmp2=*tmp;

        return tmp2;
    }
    ComplexVariable& solver::operator+(ComplexVariable& other ,complex<double> c)
    {
        ComplexVariable *tmp = new ComplexVariable;

        tmp->n1=other.n1;
        tmp->n2=other.n2;
        tmp->n3=other.n3 + c;

        ComplexVariable &tmp2=*tmp;

        return tmp2;
    }
    ComplexVariable& solver::operator+(complex<double> c, ComplexVariable& other)
    {
        ComplexVariable *tmp = new ComplexVariable;

        tmp->n1=other.n1;
        tmp->n2=other.n2;
        tmp->n3=c + other.n3;

        ComplexVariable &tmp2=*tmp;

        return tmp2;
    }

    



    ComplexVariable& ComplexVariable::operator / (ComplexVariable &other)
    {
          
            ComplexVariable *tmp = new ComplexVariable;

            tmp->n1 = this->n1 / other.n1;
            tmp->n2 = this->n2 / other.n2;
            tmp->n3 = this->n3 / other.n3;

            ComplexVariable &tmp2=*tmp;

            return tmp2;
        
    }
    ComplexVariable& solver::operator / (double n, ComplexVariable &other)
    {
        ComplexVariable *tmp = new ComplexVariable;

        tmp->n1 = n / other.n1;
        tmp->n2 = n / other.n2;
        tmp->n3 = n / other.n3;

        ComplexVariable &tmp2=*tmp;

        return tmp2;    
    }
    ComplexVariable& solver::operator / (ComplexVariable &other, double n)
    {
        if(n==0)
            throw std::invalid_argument("exception"); 

        ComplexVariable *tmp = new ComplexVariable;

        tmp->n1 = other.n1 / n;
        tmp->n2 = other.n2 / n;
        tmp->n3 = other.n3 / n;

        ComplexVariable &tmp2=*tmp;

        return tmp2;    
        
    }
    ComplexVariable& solver::operator / (complex<double> c, ComplexVariable &other)
    { 
            ComplexVariable *tmp = new ComplexVariable;

            tmp->n1 = c / other.n1;
            tmp->n2 = c / other.n2;
            tmp->n3 = c / other.n3;

            ComplexVariable &tmp2=*tmp;

            return tmp2;
        
    }
    ComplexVariable& solver::operator / (ComplexVariable &other, complex<double> c)
    {
        if(c.imag() != 0. && c.real() != 0.)
        {  
            ComplexVariable *tmp = new ComplexVariable;

            tmp->n1 = other.n1 / c;
            tmp->n2 = other.n2 / c;
            tmp->n3 = other.n3 / c;

            ComplexVariable &tmp2=*tmp;

            return tmp2;
        }
        else
        {
            throw std::invalid_argument("exception");
        }
        return other;

    }

   






