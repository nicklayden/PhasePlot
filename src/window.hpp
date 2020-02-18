#pragma once
#include <GL/glut.h>
#include <SFML/Graphics.hpp>
#include <math.h>
#include <vector>
#include <boost/program_options.hpp>
/*
	OpenGL / SFML initialization and drawing functions.


*/

void glSetup(boost::program_options::variables_map& vm, int argc, char** argv);
void cpDraw_BoundingSphere();
void cpDraw_Axes();
int PollingCenter(sf::Event event, bool& running, std::map<std::string, double>& controls, std::map<std::string, bool>& rotation);
void transformations(std::map<std::string, double>& controls, std::map<std::string, bool>& rotation);
template <class T>
void plotSurface(T (*f)(T,T),std::vector<T> x,std::vector<T> y, float r = 0.0, float g = 0.0, float b = 1.0);

template <class T>
void plotSurface_fixed(std::vector<T> x,std::vector<T> y, std::vector<std::vector<T> > z, float r=0., float =0., float b=1.0);

template<class T>
std::vector<std::vector<T> > surface(T (*f)(T,T), std::vector<T> xs, std::vector<T> ys); 



template<class T>
std::vector<std::vector<T> > surface(T (*f)(T,T), std::vector<T> xs, std::vector<T> ys) {
    // Evaluating the function on the mesh here.
    std::vector<std::vector<T> > zs;
    std::vector<T> slice;
    for (size_t i = 0; i < xs.size(); i++)
    {
        for (size_t j = 0; j < ys.size(); j++)
        {
            slice.push_back(f(xs[i],ys[j]));
        }
        zs.push_back(slice);
        slice.clear();
    }
    return zs;

}

template <class T>
void plotSurface_fixed(std::vector<T> x,std::vector<T> y, std::vector<std::vector<T> > z, float r, float g, float b) {
    // Plots the surface through sections of f(x,y)
    // This function does not call f at each step, clearly better than the other one.
    // std::vector<T> zs;
    for (size_t i = 0; i < x.size(); i++)
    {
        glBegin(GL_LINE_STRIP);
            glColor3f(r,g,b);
            for (size_t j = 0; j < y.size(); j++)
            {  
                glVertex3f(x[i],y[j],z[i][j]/250.);   
            }
        glEnd();
    }

    for (size_t i = 0; i < y.size(); i++)
    {
        glBegin(GL_LINE_STRIP);
            glColor3f(r,g,b);
            for (size_t j = 0; j < x.size(); j++)
            {   
                glVertex3f(x[j],y[i],z[j][i]/250.);   
            }
        glEnd();
    }
}

template <class T>
void plotSurface(T (*f)(T,T),std::vector<T> x,std::vector<T> y, float r, float g, float b) {
    // Plots the surface through sections of f(x,y)
    // Optimization -> Take function evaluations out of each drawing step.
    std::vector<T> z;
    for (size_t i = 0; i < x.size(); i++)
    {
        glBegin(GL_LINE_STRIP);
            glColor3f(r,g,b);
            for (size_t j = 0; j < y.size(); j++)
            {  
                if(j > 0 && z[j]*z[j-1]<0){
                    // glColor3f(1,0,0);
                    glColor3f(r,g,b);
                } else {
                    glColor3f(r,g,b);
                } 
                z.push_back(f(x[i],y[j]));
                glVertex3f(x[i]-M_PI/4 +0.1,y[j]-0.9,z[j]/500.);   
            }
            z.clear();
        glEnd();
    }

    for (size_t i = 0; i < y.size(); i++)
    {
        glBegin(GL_LINE_STRIP);
            glColor3f(r,g,b);
            for (size_t j = 0; j < x.size(); j++)
            {   
                if(j > 0 && z[j]*z[j-1]<0){
                    // glColor3f(1,0,0);
                    glColor3f(r,g,b);
                } else {
                    glColor3f(r,g,b);
                }
                z.push_back(f(x[j],y[i]));
                glVertex3f(x[j]-M_PI/4 +0.1,y[i]-0.9,z[j]/500.);   
            }
            z.clear();
        glEnd();
    }
}