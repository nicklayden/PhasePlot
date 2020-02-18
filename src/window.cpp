#include "window.hpp"



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