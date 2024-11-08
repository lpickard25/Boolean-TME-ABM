// Thie code contained herein is used for visualizing the ABM, generating and saving images and gif of the model while it runs.
// Created by Rebecca Bekker on 8/14/24.
//
#include "Environment.h"

#include <Cell.h>
#include <GL/gl.h>
#include <GLFW/glfw3.h>
#include <stb_image_write.h>
#include <cmath>
#include <vector>


/** This function creates an openGL window, and returns the pointer of that window.
 * @param width, height: these parameters define the size of the window. These are currently relative to the screen used,
 * so figures created on different machines may be different sizes.
 * @param title: Sets the title of the window.
 * @param active: this boolean toggles the visualization of the model on and off. Helpful for running on hardware
 * which do not have screens (like most HPCs).
 */

GLFWwindow* Environment::createWindow(int width, int height, const char* title, bool active) {

        std::cout<<"In create window 1"<< std::endl;
        if (active) {
                if(!glfwInit()) { // check that the library has been correctly initialized
                        std::cout<<"In create window 2"<< std::endl;
                        std::cout<<"GFLW not initialized"<<std::endl;
                        return nullptr;
                }

                GLFWwindow* win = glfwCreateWindow(width, height, title, nullptr, nullptr);

                int windowWidth, windowHeight;
                glfwGetWindowSize(win,&windowWidth,&windowHeight);

                const GLFWvidmode* vidmode = glfwGetVideoMode(glfwGetPrimaryMonitor());
                int widthPix = windowWidth;
                int heightPix = windowHeight;

                glfwSetWindowPos(win,(vidmode->width - widthPix)/2,(vidmode->height - heightPix)/2);
                glfwMakeContextCurrent(win);

                glClearColor(256,256,256,0);
                glClear(GL_COLOR_BUFFER_BIT);

                glfwDefaultWindowHints();
                glfwShowWindow(win);

                glMatrixMode(GL_PROJECTION);
                glLoadIdentity();
                glOrtho(-widthPix/2.0, widthPix/2.0, -heightPix/2.0, heightPix/2.0, -1, 1);
                glMatrixMode(GL_MODELVIEW);

                if (!win) {
                        std::cout << "Error: window not created"<<std::endl;
                        glfwTerminate();
                        return nullptr;
                }
              return win;
        }
        std::cout<<"In create window: failed completely"<< std::endl;
        return nullptr;
}

/** This function handles the drawing of the ABM.
 * @param active: this boolean toggles the visualization of the model on and off.
 * @param win: the GLFWwindow pointer returned by the createWindow function, which is the window on which drawModel draws,
 * @param scale: used to scale the visualization
 *
 * This function calls the createSubLists, drawCircle and drawCircleOutlines functions, amongst others.
 */

void Environment::drawModel(bool active, GLFWwindow* win, int scale) {
        // If visualization is switched on
        if (active) {
                glClearColor(256,256,256,0); // Set the background color to white.
                glClear(GL_COLOR_BUFFER_BIT);
                glfwMakeContextCurrent(win);

                // Sorts the cells according to the cell type (NOT STATE) into separate cell lists, so that we can layer each cellular population.
                auto sublists = createSubLists();

                // This is the order in which cells are drawn. Since Cancer cells are the largest, they're drawn first,
                // followed by M0, M1, M2, CD4 Th, Treg, CD8 active and CD8 suppressed.

                std::vector<int> drawOrder = {3,0,1,2,4,5,6,7};
                glEnable(GL_BLEND);
                glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

                for(int i : drawOrder) {
                        // Draws the cell
                        for(const auto& c :sublists[i]) {
                                // Draws and shades the circle of the cell, according to the color scheme in, outlines the cell in black
                                drawCircle(c.x[0],c.x[1],c.radius,  setColor(i), scale);
                        }
                }
                // This updates the drawing. Without this line you won't see any changes on the window.
                updateWindow(win);
        }
}


std::array<std::vector<Cell>,8> Environment::createSubLists() {
        std::array<std::vector<Cell>,8> sublists;

        for(const auto& c: cell_list) {
                sublists[c.state].push_back(c);
        }
        return sublists;

}

void Environment::updateWindow(GLFWwindow* win) {
        if(win) {
                glfwSwapBuffers(win);
                glfwPollEvents();
        }
}


std::vector<int> Environment::setColor(int cellState) {

        // Color scheme can be found at: https://personal.sron.nl/~pault/ (Figure 4: muted color scheme),
        // and is already color-blind friendly.

        /** CELL STATES
         * -1 - dead
         * 0 - M0
         * 1 - M1
         * 2 - M2
         * 3 - alive (cancer)
         * 4 - Th
         * 5 - Treg
         * 6 - active (CD8)
         * 7 - suppressed (CD8)
         */
        std::vector<int> INDIGO = {51, 34, 136};
        std::vector<int> CYAN = {136, 204, 238};
        std::vector<int> TEAL = {64, 170, 153};
        std::vector<int> GREEN = {17, 119, 51};
        std::vector<int> OLIVE = {152, 153, 51};
        std::vector<int> SAND = {221, 204, 119};
        std::vector<int> ROSE = {204, 102, 119};
        std::vector<int> WINE = {136, 34, 85};
        std::vector<int> PURPLE = {170, 68, 153};
        std::vector<int> PALE_GREY = {221, 221, 221};
        std::vector<int> BLACK = {0, 0, 0};
        std::vector<int> WHITE = {255, 255, 255};


        switch(cellState) {
                case -1: // Dead cells are shaded olive green
                        return OLIVE;
                        break;
                case 0: // M0 cells
                        return INDIGO;
                        break;
                case 1: // M1 cells
                        return CYAN;
                        break;
                case 2: // M2 cells
                        return TEAL;
                        break;
                case 3: // Cancer cells
                        return GREEN;
                        break;
                case 4: // CD4 Th cells
                        return ROSE;
                        break;
                case 5: // CD4 Foxp3 Treg cells
                        return WINE;
                        break;
                case 6: // Active CD8+ T cells
                        return PURPLE;
                        break;
                case 7: // Suppressed CD8+ T cells
                        return SAND;
                        break;
                default:
                        return WHITE;
        }
}


/** This function saves a GLFWwindow as an image
 *
 * @param filepath      location of the saved image
 * @param win           pointer to the GLFWwindow that is saved
 *
 */
void Environment::saveToJPG(const char* filepath, GLFWwindow* win) {
        // Read data from framebuffer
        int widthPix, heightPix;
        glfwGetFramebufferSize(win, &widthPix,&heightPix);
        GLsizei nrChannels = 3; // saves all three color channels
        GLsizei stride = nrChannels * widthPix;
        stride+= (stride % 4) ? (4-stride % 4) : 0;
        GLsizei bufferSize = stride * heightPix;
        std::vector<char> buffer(bufferSize);
        glPixelStoref(GL_PACK_ALIGNMENT,4);
        glReadBuffer(GL_FRONT);
        glReadPixels(0,0,widthPix,heightPix,GL_RGB, GL_UNSIGNED_BYTE,buffer.data());
        stbi_flip_vertically_on_write(true);
        if(!stbi_write_png(filepath, widthPix, heightPix, nrChannels, buffer.data(), stride)) {
                std::cout<<"Error: can't save image."<<std::endl;
        }

}


/** This function draws and shades circles which represent the cells. This is done using GL_TRIANGLE_FAN, part of the openGL library
 * @param x:            x coordinate of the center of the circle / cell
 * @param y             y coordinate of the center of the circle / cell
 * @param radius        radius of the circle / cell
 * @param COLOR         the color used for shading the circle, set by setColo
 * @param scale         scales the size of the circle
 */
void Environment::drawCircle(float x, float y, float radius, std::vector<int> COLOR, int scale) {
        // openGL draws circles by rotating triangles around the center point (x,y). A higher numTriangles makes the circle smoother, but may take longer.
        int numTriangles = 20; //

        // Draw the and shade the cell
        glColor3ub(COLOR[0],COLOR[1],COLOR[2]);
        glBegin(GL_TRIANGLE_FAN);
        glVertex2f(x/scale,y/scale); // Center of the circle

        // Calculates the vertices of the approximating triangles
         for (int ii = 0; ii <= numTriangles; ii++)   {
                 float theta = 2.0f * 3.1415926f * float(ii) / float(numTriangles);//get the current angle
                 float xNew = (radius/scale) * cosf(theta);//calculate the x component
                 float yNew = (radius/scale) * sinf(theta);//calculate the y component
                 glVertex2f(x/scale + xNew, y/scale + yNew);//output vertex
         }
        glFlush(); // Flushes the drawing buffer (draw to the window)
        glEnd(); // Ends the GL_TRIANGLE_FAN drawing

        // Draw the cell outlines
        // Note: The rationale for drawing the outline after the circle is to make sure the outline is visible, but
        // this doesn't work perfectly. I'm still trying to figure out why.
        glLineWidth(10);
        glColor3ub(0,0,0);
        glBegin(GL_TRIANGLES); //
        for (int ii = 0; ii <= numTriangles; ii++)   {
                float theta = 2.0f * 3.1415926f * float(ii) / float(numTriangles);//get the current angle
                float xNew = radius/scale * cosf(theta);//calculate the x component
                float yNew = radius/scale * sinf(theta);//calculate the y component
                glVertex2f(x/scale + xNew, y/scale + yNew);//output vertex
        }
        glFlush(); // Flushes the drawing buffer which forces all preceding GL commands to be executed. That is, everything get's drawn
        glEnd(); // Ends the drawing
}


/** This function draws the outlines of cells, similar to the functionality of drawCircle
 * @param x
 * @param y
 * @param radius
 * @param scale
 *
 * Note: This code is currently not used.
 */
void Environment::drawCellOutlines(float x, float y, float radius,int scale) {
        int numPoints = 60;
        glColor3ub(0,0,0);
        glBegin(GL_TRIANGLES);
        glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
        glEnable(GL_LINE_SMOOTH);
        glLineWidth(10);

        for (int ii = 0; ii <= numPoints; ii++)   {
                float theta = 2.0f * 3.1415926f * float(ii) / float(numPoints);//get the current angle
                float xNew = radius/scale * cosf(theta);//calculate the x component
                float yNew = radius/scale * sinf(theta);//calculate the y component
                glVertex2f(x/scale + xNew, y/scale + yNew);//output vertex
        }
        glFlush(); // Flushes the drawing buffer
        glEnd(); // Ends the drawing
        glDisable(GL_LINE_SMOOTH);
}





