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

 


/***
 * Function Declarations and Templates for testing
 * 
 * 
 */
void cpDraw_Axes();
void cpDraw_BoundingSphere();
std::vector<std::vector<double> > transpose_copy(std::vector<std::vector<double> > data);
int PollingCenter(sf::Event event, bool& running, std::map<std::string, double>& controls, std::map<std::string, bool>& rotation);
void transformations(std::map<std::string, double>& controls, std::map<std::string, bool>& rotation);
void glSetup(boost::program_options::variables_map& vm, int argc, char** argv);

template <class T>
void plotSurface(T (*f)(T,T),std::vector<T> x,std::vector<T> y, float r = 0.0, float g = 0.0, float b = 1.0);

inline double f(double z, double eta) {
    return exp(-z*z - eta*eta);
}

inline double zero_f(double r, double k){
    return 0;
}

inline double alan_f(double r, double k) {
    double a,b,c,d,e,f,g,h;
    a = sinh(k*(M_PI - r) + sinh(k*r));
    b = cosh(k*(M_PI - r) - cosh(k*r));
    c = sinh(k*(M_PI/2. + r) + sinh(k*(M_PI/2. - r)));
    d = cosh(k*(M_PI/2. + r) - cosh(k*(M_PI/2. - r)));

    e = -(1./sin(r));
    f = -cos(r)/(sin(r)*sin(r));
    g = 1./cos(r);
    h = tan(r);

    return e*(  f*a + k*b + g*( h*c + k*d )  );
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




int main(int argc, char** argv ){  

    /* #####################################################################################
        LOCAL VARIABLE DECLARATIONS AND PROGRAM SETTINGS
       #####################################################################################
    */
    // Create SDL window instance
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
    std::vector<double> xs,ys,zs;
    double begin,end,step,r0,rn,k0,kn;
    uint N;
    begin = -5;
    end = 5;
    N = 40;
    step = (end - begin)/N;

    // for (size_t i = 0; i < N; i++)
    // {
    //     xs.push_back(begin + step*i);
    //     ys.push_back(begin + step*i);
    // }
    
    r0 = M_PI/4. - M_PI/4 + 0.1;
    rn = M_PI/4. + M_PI/4;
    k0 = 1.01;
    kn = 2;

    for (size_t i = 0; i < N; i++)
    {
        xs.push_back(r0 + ((rn-r0)/N)*i );
        ys.push_back(k0 + ((kn-k0)/N)*i);
    }
    




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
        #####################################################################################
    */
    boost::program_options::variables_map vm;
    status = config_mapping(argc, argv, vm);
    std::cout << "Mapping status: " << status << std::endl;
    r = vm["init.conds.radius"].as<double>();
    dr = vm["init.conds.dr"].as<double>();
    tmax = vm["init.conds.thetamax"].as<double>();
    xcenter = vm["init.conds.xcenter"].as<double>();
    ycenter = vm["init.conds.ycenter"].as<double>();
    numinit = vm["init.conds.numinit"].as<int>();
    sysfile = vm["system.file"].as<std::string>();
    sysmethod = vm["system.method"].as<std::string>();
    maxtime = vm["system.timemax"].as<double>();
    mintime = vm["system.timemin"].as<double>();
    dt = vm["system.dt"].as<double>();
    bkg_r = vm["gui.bkg_r"].as<float>();
    bkg_g = vm["gui.bkg_g"].as<float>();
    bkg_b = vm["gui.bkg_b"].as<float>();
    bkg_alpha = vm["gui.bkg_alpha"].as<float>();
    gui_fovy = vm["gui.fovy"].as<float>();
    gui_aspect = vm["gui.aspect"].as<float>();
    gui_znear = vm["gui.znear"].as<float>();
    gui_zfar = vm["gui.zfar"].as<float>();
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
            plotSurface(alan_f,xs,ys);
            plotSurface(zero_f,xs,ys,0.,1.,0.);

            // Display everything we've done to the SFML instance
            mainwin.display();
            // This swaps the opengl buffer to sfml.

        } // end while running loop

} // end main

void glSetup(boost::program_options::variables_map& vm, int argc, char** argv) {
    /**
     * Initializes the opengl context for drawing, loads all configurations from .ini,
     * sets up the proper initial transformations
    */
    std::string sysmethod, sysfile;
    int numinit, status;
    double r, dr, tmax, xcenter, ycenter,dt,maxtime,mintime;
    float bkg_alpha, bkg_r, bkg_g, bkg_b, \
          gui_aspect, gui_zfar, gui_znear, \
          gui_fovy, gui_camera_dist,gui_camera_dist_x,gui_camera_dist_y;
    status = config_mapping(argc, argv, vm);
    std::cout << "Mapping status: " << status << std::endl;
    r = vm["init.conds.radius"].as<double>();
    dr = vm["init.conds.dr"].as<double>();
    tmax = vm["init.conds.thetamax"].as<double>();
    xcenter = vm["init.conds.xcenter"].as<double>();
    ycenter = vm["init.conds.ycenter"].as<double>();
    numinit = vm["init.conds.numinit"].as<int>();
    sysfile = vm["system.file"].as<std::string>();
    sysmethod = vm["system.method"].as<std::string>();
    maxtime = vm["system.timemax"].as<double>();
    mintime = vm["system.timemin"].as<double>();
    dt = vm["system.dt"].as<double>();
    bkg_r = vm["gui.bkg_r"].as<float>();
    bkg_g = vm["gui.bkg_g"].as<float>();
    bkg_b = vm["gui.bkg_b"].as<float>();
    bkg_alpha = vm["gui.bkg_alpha"].as<float>();
    gui_fovy = vm["gui.fovy"].as<float>();
    gui_aspect = vm["gui.aspect"].as<float>();
    gui_znear = vm["gui.znear"].as<float>();
    gui_zfar = vm["gui.zfar"].as<float>();
    gui_camera_dist = vm["gui.camera_dist"].as<float>();
    // preparing opengl surface. defining default window color.
    glClearDepth(1.f);
    glClearColor(bkg_r, bkg_g, bkg_b, bkg_alpha);
    glEnable(GL_DEPTH_TEST); // Allows depth testing for objects at various distances from camera
    glDepthMask(GL_TRUE); // Enables the depth mask for proper drawing of distance objects
    glEnable(GL_BLEND); 
    glHint(GL_LINE_SMOOTH, GL_NICEST); // Tells opengl to create smooth lines
    // setup perspective projection and the cameras position:
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(gui_fovy, gui_aspect, gui_znear, gui_zfar); // fovy, aspect ratio, znear, zfar // objects farther than zfar or closer than znear will NOT be rendered.
    glPointSize(5); // size of point sprites on the window. default is 1 pixel.
    glLineWidth(2);
    glutInit(&argc, argv);
}

void cpDraw_BoundingSphere() {
    // Draws 3 circles. One per plane in 3D, as a reference object.
    // Drawing circles around the axis for spatial reference.
    std::vector<float> xcirc,ycirc;
    float th;
    for (int i = 0; i < 80; i++) {
        th = 2.*M_PI*i/80.;
        xcirc.push_back(cos(th));
        ycirc.push_back(sin(th));
    }
    glBegin(GL_LINE_LOOP);
        for (int i = 0; i < xcirc.size(); i++) {
            glColor3f(1,0,0);
            glVertex3f(xcirc[i], ycirc[i], 0);
        }
    glEnd();
    glBegin(GL_LINE_LOOP);
        glColor3f(0,0.5,1);
        for (int i = 0; i < xcirc.size(); i++) {
            glVertex3f(xcirc[i], 0, ycirc[i]);
        }
    glEnd();
    glBegin(GL_LINE_LOOP);
        glColor3f(1,0,1);
        for (int i = 0; i < xcirc.size(); i++) {
            glVertex3f(0, xcirc[i], ycirc[i]);
        }
    glEnd();
}

void cpDraw_Axes() {
    /****************************************************************************************
     * 
     *                      AXES AND POINT DRAWING
     *  This section has no numerics. Only drawing xyz axes lines, and axes unit circles 
     *  for reference in the diagram. Positive axes are shown with colours according to:
     *      X-Axis : Red
     *      Y-Axis : Green
     *      Z-Axis : Blue
     * **************************************************************************************/


    // Draw XYZ Axes on screen.
    glBegin(GL_LINE_LOOP);
        // X AXIS in positive direction.
        glColor3f(1,0,0); // red
        glVertex3f(0,0,0);
        glVertex3f(1,0,0);
    glEnd();

    glBegin(GL_LINE_LOOP);
    // Y AXIS in positive direction.
        glColor3f(0,1,0); // green
        glVertex3f(0,0,0);
        glVertex3f(0,1,0);
    glEnd();

    glBegin(GL_LINE_LOOP);
    // Z AXIS in positive direction.
        glColor3f(0,0,1); // blue
        glVertex3f(0,0,0);
        glVertex3f(0,0,1);
    glEnd();

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
                    glColor3f(1,0,0);
                } else {
                    glColor3f(r,g,b);
                } 
                z.push_back(f(x[i],y[j]));
                glVertex3f(x[i],y[j],z[j]/250.);   
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
                    glColor3f(1,0,0);
                } else {
                    glColor3f(r,g,b);
                }
                z.push_back(f(x[j],y[i]));
                glVertex3f(x[j],y[i],z[j]/250.);   
            }
            z.clear();
        glEnd();
    }
}


int PollingCenter(sf::Event event, bool& running, std::map<std::string, double>& controls, std::map<std::string, bool>& rotation) {
        /**
         * Need to fix this commented piece here, the transformation order may not be correct.
         * It leaves a black screen on the drawing window.
        */
        // if (event.type == sf::Event::Resized) {
        //     // adjust the viewport when the window is resized, adjusting
        //     // the aspect ratio so we don't stretch the image.
        //     float new_aspect = (float)mainwin.getSize().x/(float)mainwin.getSize().y;
        //     // Adjusts the size of the viewport to reflect window resizing.
        //     glViewport(0, 0, event.size.width, event.size.height);
        //     // this is old opengl, switch to projection matrix mode,
        //     // load identity so we dont stack projections by mistake,
        //     // then change our perspective to the new aspect ratio
        //     // then finalize by switching back to modelview matrix mode
        //     // so that the drawings arent fucked up after.
        //     glMatrixMode(GL_PROJECTION);
        //     glLoadIdentity();
        //     gluPerspective((float)controls["gui_fovy"], (float)new_aspect, (float)controls["gui_znear"], (float)controls["gui_zfar"]);
        //     glMatrixMode(GL_MODELVIEW);
        // }
        if ((event.type == sf::Event::KeyPressed) && (event.key.code == sf::Keyboard::Right)){
            controls["rotspeed"] += 1;
        }
        else if ((event.type == sf::Event::KeyPressed) && (event.key.code == sf::Keyboard::Left)){
            controls["rotspeed"] -= 1;
        }
        // X ROTATIONS
        else if ((event.type == sf::Event::KeyPressed) && (event.key.code == sf::Keyboard::Numpad7)){
            rotation["rotatex"]=!rotation["rotatex"]; // toggle continuous rotation in x
        }
        else if ((event.type == sf::Event::KeyPressed) && (sf::Keyboard::isKeyPressed(sf::Keyboard::J))){
            controls["anglex"] += 0.1; // hold down positive rotation in x
        }
        else if ((event.type == sf::Event::KeyPressed) && (sf::Keyboard::isKeyPressed(sf::Keyboard::U))){
            controls["anglex"] -= 0.1; // hold down positive rotation in x
        }
        // Y ROTATIONS
        else if ((event.type == sf::Event::KeyPressed) && (event.key.code == sf::Keyboard::Numpad8)){
            rotation["rotatey"]=!rotation["rotatey"]; // toggle continuous rotation in y
        }
        else if ((event.type == sf::Event::KeyPressed) && (sf::Keyboard::isKeyPressed(sf::Keyboard::I))){
            controls["angley"] += 0.1; // hold down positive rotation in y
        }
        else if ((event.type == sf::Event::KeyPressed) && (sf::Keyboard::isKeyPressed(sf::Keyboard::K))){
            controls["angley"] -= 0.1; // hold down positive rotation in y
        }
        // Z ROTATIONS
        else if ((event.type == sf::Event::KeyPressed) && (event.key.code == sf::Keyboard::Numpad9)){
            rotation["rotatez"]=!rotation["rotatez"]; // toggle continuous rotation in z
        }
        else if ((event.type == sf::Event::KeyPressed) && (sf::Keyboard::isKeyPressed(sf::Keyboard::O))){
            controls["anglez"] += 0.1; // hold down positive rotation in z
        }
        else if ((event.type == sf::Event::KeyPressed) && (sf::Keyboard::isKeyPressed(sf::Keyboard::L))){
            controls["anglez"] -= 0.1; // hold down negative rotation in z
        } 
        else if ((event.type == sf::Event::KeyPressed) && (event.key.code == sf::Keyboard::Up)){
            controls["gui_camera_dist"] -= 0.2; // z distance to camera
        }
        else if ((event.type == sf::Event::KeyPressed) && (event.key.code == sf::Keyboard::Down)){
            controls["gui_camera_dist"] += 0.2; // z distance to camera
        }
        else if ((event.type == sf::Event::KeyPressed) && (event.key.code == sf::Keyboard::PageUp)){
            controls["gui_camera_dist_x"] -= 0.2; // x distance to camera
        }
        else if ((event.type == sf::Event::KeyPressed) && (event.key.code == sf::Keyboard::PageDown)){
            controls["gui_camera_dist_x"] += 0.2; // x distance to camera
        }
        else if ((event.type == sf::Event::KeyPressed) && (event.key.code == sf::Keyboard::Home)){
            controls["gui_camera_dist_y"] -= 0.2; // y distance to camera
        }
        else if ((event.type == sf::Event::KeyPressed) && (event.key.code == sf::Keyboard::End)){
            controls["gui_camera_dist_y"] += 0.2; // y distance to camera
        }

    return 0;
}


void transformations(std::map<std::string, double>& controls, std::map<std::string, bool>& rotation) {
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(-controls["gui_camera_dist_x"],-controls["gui_camera_dist_y"],-controls["gui_camera_dist"]);
    glRotatef(-75.,1.,0.,0.); // this rotates to put XY plane into perspective
    if (controls["rotatex"]) {
        controls["anglex"] += .05;
    }
    if (controls["rotatey"]) {
        controls["angley"] += .05;
    }
    if (controls["rotatez"]) {
        controls["anglez"] += .05;
    }
    // opengl rotation functions.
    glRotatef(controls["anglex"]*controls["rotspeed"],1.0f,0.f,0.f);
    glRotatef(controls["angley"]*controls["rotspeed"],0.0f,1.0f,0.0f);
    glRotatef(controls["anglez"]*controls["rotspeed"],0.0f,0.0f,1.0f);
}