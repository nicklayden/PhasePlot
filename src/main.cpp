/*
    PhasePlot Version 1.0
    Author: Nick Layden, 2019 (October)
        Dependencies:
            OpenGL: for drawing and displaying. (freeglut, glut, glu, glew, gl)
            SFML  : for handling the draw window and keyboard inputs for motion.
            Boost : for Program configuration files, ODE integrating etc. 

    There are a considerable number of optimizations to be done, so just compile with -O3
    
*/
// Library requirements
// OpenGL/glut
#include <GL/glut.h>
// Standard C++
#include <vector>
#include <math.h>
#include <random>
#include <iostream>
#include <map>

// SFML/ OpenGL  - Used for keyboard event handling and opengl configuration
#include <SFML/Graphics.hpp>

#include "cmdline.hpp"
#include "window.hpp"


inline double f(double z, double eta) {
    return exp(-z*z - eta*eta);
}

inline double zero_f(double r, double k){
    return 0;
}


inline double alan_f(double r, double k) {
    // Alans specific prob. Just wanted a graph showing the function
    // crossing f(r,k)=0 for r near pi/4 and k = (1,2)
    double a,b,c,d,e,f,g,h;
    a = sinh(k*(M_PI - r)) + sinh(k*r);
    b = cosh(k*(M_PI - r)) - cosh(k*r);
    c = sinh(k*(M_PI/2. + r)) + sinh(k*(M_PI/2. - r));
    d = cosh(k*(M_PI/2. + r)) - cosh(k*(M_PI/2. - r));

    e = -(1./sin(r));
    f = -cos(r)/(sin(r)*sin(r));
    g = 1./cos(r);
    h = tan(r);

    return e*(  f*a + k*b) + g*( h*c + k*d )  ;
}

inline double alan_f_A(double r, double k) {
    // Alans function
    // THIS IS A(r) 
    double a,b,c,d,e,f;
    a = sinh(k*M_PI);
    b = cosh(k*r);
    c = cosh(k*M_PI);
    d = sinh(k*r);
    e = -cos(r)/(sin(r)*sin(r));
    f = k/sin(r);

    return e*(a*b - c*d + d) + f*(a*d - c*b + b);
}

inline double alan_f_Gr(double r, double k) {
    // Alan's function G(r),r
    // Depends on alan_f_A for calculating A(r2)
    double a,b,c,d,r2;

    a = sqrt(2)*cos(r);
    b = sqrt(1. - 0.5*sin(r)*sin(r));
    c = a/b;
    r2 = acos(sin(r)*(-1/sqrt(2)));

    d = alan_f_A(r2,k);
    return c*d;
}

inline double alan_theta(double r, double k) {
    // Alan's function sqrt(theta) = F,r + G,r
    // Depends on alan_f_Gr, alan_f
    return alan_f_A(r,k) + alan_f_Gr(r,k);
}

inline double alan_ar_test(double r, double k) {
    return alan_f_A(r,k) - alan_f_A(M_PI/2. - r, k); 
}


inline double plane(double x, double y) {
    //eqn for a plane: ax + by + cz = d
    // as f(x,y) -> z = (d - ax - by)/c
    double a,b,c,d;
    a = 1.0;
    b = 3.0;
    c = 2.0;
    d = 1.0;

    return (d - a*x - b*y)/c;
}

inline double paraboloid(double x, double y) {
    return x*x + y*y;
}


class func {
    public: 
        func(double xi):t(xi) {};
        double t,a,b;
        double operator()(double x, double y) {
            a = sin(t);
            b = cos(t);
            return exp(-(x-b)*(x-b) - (y-a)*(y-a));
        }
};


class alan_func {
    public: 
        alan_func(double xi):t(xi) {};
        double t,a,b;
        double operator()(double x, double y) {
            a = sin(t);
            b = cos(t);
            return exp(-(x-b)*(x-b) - (y-a)*(y-a));
        }
};


std::vector<double> uniform_grid(double a, double b, size_t N) {
    // Create a uniform 1-D grid [a,b]
    std::vector<double> grid;
    for (int i = 0; i < N; ++i)
    {
        grid.push_back(a + i*(b-a)/N);
    }
    return grid;
}

std::vector<double> nonuniform_grid(double a, double b, size_t N, double (*param)(double,double)) {
    // Create a non-uniform grid on [a,b], using a parameterization.
}


inline double alan_const_k(double r, double k) {
    return alan_f_A(r,1.0);
}



int main(int argc, char** argv ){  

    /* #####################################################################################
        LOCAL VARIABLE DECLARATIONS AND PROGRAM SETTINGS
       #####################################################################################
    */
    // Create SFML window instance
    sf::ContextSettings settings;
    settings.depthBits = 24; settings.majorVersion = 3; settings.minorVersion = 0;
    sf::RenderWindow mainwin(sf::VideoMode(800,800), "Phase Space", sf::Style::Default, settings);

    /*************************************
     *  test plotting a custom surface.
     *  need to define x and y meshgrid, and a function to plot.
     *  in the form z=f(x,y)
     * 
     * Idea: Draw each line (slice) of the function per loop
     *       to fill up the space.
     * ***********************************/
    std::vector<std::vector<double> > zs;
    std::vector<double> xs,ys,slice;
    double begin,end,step,r0,rn,k0,kn;
    uint N=100;
    

    // Alan's 2 Variable system domain
    r0 = -3.;
    rn = 3.;
    k0 = 0.9;
    kn = 1.1;

    // Mesh grid + function on the grid to be plotted in 3D.
    xs = uniform_grid(r0,rn,N);
    ys = uniform_grid(k0,kn,N);
    zs = surface(alan_ar_test,xs,ys);
    
    // // Evaluating the function on the mesh here.
    // for (size_t i = 0; i < xs.size(); i++)
    // {
    //     for (size_t j = 0; j < ys.size(); j++)
    //     {
    //         slice.push_back(alan_f_Gr(xs[i],ys[j]));
    //     }
    //     zs.push_back(slice);
    //     slice.clear();
    // }
        
    // for (int i = 0; i < N; ++i)
    // {
    //     std::cout << "Bracket: r= " << r0 + i*(rn-r0)/N << "  f(r,k)= " << alan_ar_test(r0 + i*(rn-r0)/N, 1.2) << std::endl;
    // }



    // Camera rotation flags and properties.
    bool rotate = false;
    bool rotatex = false; 
    bool rotatey = false;
    bool rotatez = false;
    float angle = 0;
    float anglex = 90;
    float angley = 0;
    float anglez = 0;
    float rotspeed = 10;
    double lambda; // Scalar field parameter.



    // Initial conditions, dynamical system settings, and gui settings for the window.
    int numinit, status;
    double r, dr, tmax, xcenter, ycenter,dt,maxtime,mintime;
    float bkg_alpha, bkg_r, bkg_g, bkg_b, \
          gui_aspect, gui_zfar, gui_znear, \
          gui_fovy, gui_camera_dist,gui_camera_dist_x,gui_camera_dist_y;
    std::string sysmethod, sysfile;
    gui_camera_dist_x = 0.f;
    gui_camera_dist_y = 0.f;

    // Configuration configs(argc,argv,std::string("config.ini"));

    std::map<std::string, bool> Rotation_flags;
    Rotation_flags["rotate"] = rotate;
    Rotation_flags["rotatex"] = rotatex;
    Rotation_flags["rotatey"] = rotatey;
    Rotation_flags["rotatez"] = rotatez;


    
    std::map<std::string, double> Camera_controls;
    Camera_controls["angle"] = angle;
    Camera_controls["anglex"] = anglex;
    Camera_controls["angley"] = angley;
    Camera_controls["anglez"] = anglez;
    Camera_controls["rotspeed"] = rotspeed;
    Camera_controls["gui_camera_dist_x"] = gui_camera_dist_x;
    Camera_controls["gui_camera_dist_y"] = gui_camera_dist_y;
    /* ######################################################################################
        LOAD PROGRAM CONFIGURATION FROM CONFIG.INI!!!!
    
        need to capture variables from the config file, and assign them in this block of main.
        NEED TO REMOVE THIS SHIT EVENTUALLY
        #####################################################################################
    */

   /*
    This block needs to be removed from main!
   */
    boost::program_options::variables_map vm;
    status = config_mapping(argc, argv, vm);
    gui_fovy = vm["gui.fovy"].as<float>();
    gui_camera_dist = vm["gui.camera_dist"].as<float>();
    Camera_controls["gui_camera_dist"] = gui_camera_dist;
    Camera_controls["gui_fovy"] = gui_fovy;

    glSetup(vm,argc,argv);


    /**
     * Event handling on the main drawing window for display.
    */
    bool running = true; // This is used so that the loop can end and opengl frees resources properly!!
    while (running) {
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            // Handle events through the SFML interface for the window. (keyboard press, closing, etc.)
            sf::Event event;
            // PollingCenter(&mainwin, running, Camera_controls,Rotation_flags);
            while (mainwin.pollEvent(event)) {
                if (event.type == sf::Event::Closed){
                    running = false;
                    break;
                } 
                status = PollingCenter(event, running, Camera_controls, Rotation_flags);
            }
            transformations(Camera_controls,Rotation_flags);
            // Draw the axes for a reference grid.
            cpDraw_Axes();
            // Plot the things we want to see
            plotSurface(alan_const_k,xs,ys);
            plotSurface(zero_f,xs,ys,0.5,0.,0.5);
            // plotSurface(paraboloid,xs,ys,1.0,0.,0.2);
            // plotSurface_fixed(xs,ys,zs);

            // Display everything we've done to the SFML instance
            mainwin.display();
            // This swaps the opengl buffer to sfml.

        } // end while running loop

} // end main
