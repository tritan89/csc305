// C++ include
#include <iostream>
#include <string>
#include <vector>

// Utilities for the Assignment
#include "utils.h"

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"

// Shortcut to avoid Eigen:: everywhere, DO NOT USE IN .h
using namespace Eigen;

float quadSolver(const Vector3d &d, const Vector3d &e, const Vector3d &c, float r ){
    const Vector3d b = e - c;
    float discrim = d.dot(b)*d.dot(b) - d.dot(d)*(b.dot(b)- r*r);
    float root1 , root2;
    if (discrim < 0 ) return false;
    else if (discrim == 0 ) return  (-d.dot(b)) /(d.dot(d));
    else {
        root1 = (-d.dot(b) - sqrt(discrim))/(d.dot(d));
        root2 = (-d.dot(b) + sqrt(discrim))/(d.dot(d));
        if (root1 < 0 && root2 < 0)
           return false;

     if (root1 < root2)
            return root1;
        else return root2;
}


    
}

void raytrace_sphere()
{
    std::cout << "Simple ray tracer, one sphere with orthographic projection" << std::endl;

    const std::string filename("sphere_orthographic1.png");
    MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);



    // The camera is orthographic, pointing in the direction -z and covering the
    // unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1); //C
    const Vector3d x_displacement(2.0 / C.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / C.rows(), 0);

    const Vector3d sphere_center(0, -10, -3);
    const double sphere_radius = 0.9;
    // Single light source
    const Vector3d light_position(-1, 1, 1);

    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = sphere_center + double(i) * x_displacement + double(j) * y_displacement;

            // Prepare the ray
            const Vector3d ray_origin = camera_origin; //O
            const Vector3d ray_direction = pixel_center- camera_origin; //D





            float inter = quadSolver(ray_direction, ray_origin, sphere_center, sphere_radius);
            bool hit;
            if ( inter == false){
                hit = false;
            }
            else{
                hit = true;
            }

            if (hit)
            {
                // The ray hit the sphere, compute the exact intersection point
                Vector3d ray_intersection = ray_origin + ray_direction*inter;

                // Compute normal at the intersection point
                Vector3d ray_normal = ((ray_intersection - sphere_center)/sphere_radius).normalized();

                // Simple diffuse model
                C(i, j) = (light_position - ray_intersection).normalized().transpose() * ray_normal;

                // Clamp to zero
                C(i, j) = std::max(C(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(C, C, C, A, filename);
}

void raytrace_perspective()
{
  std::cout << "Simple ray tracer, one parallelogram with perspective projection" << std::endl;

    const std::string filename("plane_perspective.png");
    MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is perspective, pointing in the direction -z and covering the unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / C.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / C.rows(), 0);

    // TODO: Parameters of the parallelogram (position of the lower-left corner + two sides)
    const Vector3d pgram_origin(-0.5, -0.5, 0);
    const Vector3d pgram_u(0, 0.7, -10);
    const Vector3d pgram_v(1, 0.4, 0);


    


    // Single light source
    const Vector3d light_position(-1, 1, 1);

    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // Prepare the ray
            const Vector3d ray_origin = camera_origin;
            const Vector3d ray_direction = pixel_center- camera_origin;
            
            //variables for intersection 
            double a = pgram_u[0];
            double b = pgram_u[1];
            double c = pgram_u[2];

            double d = pgram_v[0];
            double e = pgram_v[1];
            double f = pgram_v[2];

            double g = -ray_direction[0];
            double h = -ray_direction[1];
            double i_ = -ray_direction[2];

            double j_ = -(pgram_origin[0] - ray_origin[0]);
            double k = -(pgram_origin[1] - ray_origin[1]);
            double l = -(pgram_origin[2] - ray_origin[2]);

            // // solve for each varible B Y T 


            double M  = a*(e*i_ -h*f)+ b*(g*f-d*i_) + c*(d*h-e*g);
            double theta =-( f* (a*k -j_*b)+ e*(j_*c-a*l) + d*(b*l-k*c))/M;
            double lam =( i_*(a*k-j_*b) + h*(j_*c-a*l) + g*(b*l-k*c))/M;
            double beta  = (j_*(e*i_-h*f) +k*(g*f-d*i_)+ l*(d*h-e*g))/M;
           
           
            bool intersect = false;
            if (theta >=0 && (lam<= 1 && lam>=0) && (beta<= 1 && beta>=0)){
                intersect = true;
            }




            // TODO: Check if the ray intersects with the parallelogram
            
            if (intersect)
            {
               
                
                Vector3d ray_intersection= ray_origin + theta*ray_direction;

                // TODO: Compute normal at the intersection point
                Vector3d ray_normal = pgram_v.normalized().cross(pgram_u.normalized());

                // Simple diffuse model
                C(i, j) = (light_position - ray_intersection).normalized().transpose() * ray_normal;

                // Clamp to zero
                C(i, j) = std::max(C(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(C, C, C, A, filename);
}

void  raytrace_parallelogram()
{
    
    std::cout << "Simple ray tracer, one parallelogram with orthographic projection" << std::endl;

    const std::string filename("plane_orthographic.png");
    MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is perspective, pointing in the direction -z and covering the unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / C.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / C.rows(), 0);

    // TODO: Parameters of the parallelogram (position of the lower-left corner + two sides)
    const Vector3d pgram_origin(-0.5, -0.5, 0);
    const Vector3d pgram_u(0, 0.7, -10);
    const Vector3d pgram_v(1, 0.4, 0);


    


    // Single light source
    const Vector3d light_position(-1, 1, 1);

    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // Prepare the ray
            const Vector3d ray_origin = pixel_center;
            const Vector3d ray_direction = camera_view_direction;
            
            //variables for intersection 
            double a = pgram_u[0];
            double b = pgram_u[1];
            double c = pgram_u[2];

            double d = pgram_v[0];
            double e = pgram_v[1];
            double f = pgram_v[2];

            double g = -ray_direction[0];
            double h = -ray_direction[1];
            double i_ = -ray_direction[2];

            double j_ = -(pgram_origin[0] - ray_origin[0]);
            double k = -(pgram_origin[1] - ray_origin[1]);
            double l = -(pgram_origin[2] - ray_origin[2]);

            // // solve for each varible B Y T 


            double M  = a*(e*i_ -h*f)+ b*(g*f-d*i_) + c*(d*h-e*g);
            double theta =-( f* (a*k -j_*b)+ e*(j_*c-a*l) + d*(b*l-k*c))/M;
            double lam =( i_*(a*k-j_*b) + h*(j_*c-a*l) + g*(b*l-k*c))/M;
            double beta  = (j_*(e*i_-h*f) +k*(g*f-d*i_)+ l*(d*h-e*g))/M;
           
           
            bool intersect = false;
            if (theta >=0 && (lam<= 1 && lam>=0) && (beta<= 1 && beta>=0)){
                intersect = true;
            }




            // TODO: Check if the ray intersects with the parallelogram
            
            if (intersect)
            {
               
                
                Vector3d ray_intersection= ray_origin + theta*ray_direction;

                // TODO: Compute normal at the intersection point
                Vector3d ray_normal = pgram_v.normalized().cross(pgram_u.normalized());

                // Simple diffuse model
                C(i, j) = (light_position - ray_intersection).normalized().transpose() * ray_normal;

                // Clamp to zero
                C(i, j) = std::max(C(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(C, C, C, A, filename);
}

void raytrace_shading()
{
    std::cout << "Simple ray tracer, one sphere with different shading" << std::endl;

    const std::string filename("shading.png");
    MatrixXd R = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd G = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd B = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is perspective, pointing in the direction -z and covering the unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / A.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / A.rows(), 0);

    //Sphere setup
    const Vector3d sphere_center(0, 0, 0);
    const double sphere_radius = 0.9;

    //material params
    const Vector3d diffuse_color(1, 0, 1);
    const double specular_exponent = 100;
    const Vector3d specular_color(0., 0, 1);

    // Single light source
    const Vector3d light_position(-1, 1, 1);
    double ambient = 0.1;

    for (unsigned i = 0; i < R.cols(); ++i)
    {
        for (unsigned j = 0; j < R.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // Prepare the ray
            const Vector3d ray_origin = pixel_center;
            const Vector3d ray_direction = camera_view_direction;

            // Intersect with the sphere
            // NOTE: this is a special case of a sphere centered in the origin and for orthographic rays aligned with the z axis
            Vector2d ray_on_xy(ray_origin(0), ray_origin(1));
            const double sphere_radius = 0.9;

            if (ray_on_xy.norm() < sphere_radius)
            {
                // The ray hit the sphere, compute the exact intersection point
                Vector3d ray_intersection(
                    ray_on_xy(0), ray_on_xy(1),
                    sqrt(sphere_radius * sphere_radius - ray_on_xy.squaredNorm()));

                // Compute normal at the intersection point
                Vector3d ray_normal = ray_intersection.normalized();

                // TODO: Add shading parameter here

                // need l a normal vector to light position
                // need v a unit vector to view position
                // need n surface normal vector
                Vector3d l =  (light_position - ray_intersection).normalized();
                Vector3d v =   (camera_origin - ray_intersection).normalized();
                Vector3d h = ((l+v)/(l+v).norm()).normalized();
                const double diffuse = std::max(0.0, (light_position - ray_intersection).normalized().dot(ray_normal));
                const double specular = std::max(0.0, ray_normal.dot(h));
                const double phong = pow(specular, specular_exponent);
                // Simple diffuse model
                R(i, j) = ambient + diffuse*diffuse_color(0) + phong*specular_color(0);
                // Clamp to zero
                R(i, j) = std::max(R(i, j), 0.);


                G(i, j) = ambient + diffuse*diffuse_color(1) + phong*specular_color(1);

                // Clamp to zero
                G(i, j) = std::max(G(i, j), 0.);

                B(i, j) = ambient + diffuse*diffuse_color(2) + phong*specular_color(2);

                // Clamp to zero
                B(i, j) = std::max(B(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(R, G, B, A, filename);
}

int main()
{
    raytrace_sphere();
    //raytrace_parallelogram();
    //raytrace_perspective();
    raytrace_shading();

    return 0;
}
