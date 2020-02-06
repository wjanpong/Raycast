
#include "include/Angel.h"
#include <stdio.h>
#include <GL/glut.h>
#include <math.h>
#include "global.h"
#include "sphere.h"
#include <cstdlib>
#include <vector>
// Global variables
//
using namespace std;

extern int win_width;
extern int win_height;

extern GLfloat frame[WIN_HEIGHT][WIN_WIDTH][3];  

extern float image_width;
extern float image_height;

extern Point eye_pos;
extern float image_plane;
extern RGB_float background_clr;
extern RGB_float null_clr;

extern Spheres *scene;

// light 1 position and color
extern Point light1;
extern float light1_intensity[3];

// global ambient term
extern float global_ambient[3];

// light decay parameters
extern float decay_a;
extern float decay_b;
extern float decay_c;

extern int shadow_on; //extern in raycast.cpp
extern int step_max;
extern int shadow_on;
extern int reflection_on;
extern int refraction_on;
extern int chessboard_on;
extern int stochastic_on;
extern int super_sampling_on;

float board_ambient[3]={0.4,0.4,0.4};
float board_specular[3]={1.0,1.0,1.0};
float board_shiness = 10;
int max_rays = MAX_RAYS;
int generated_ray = 0;

extern float refractive_index;

/////////////////////////////////////////////////////////////////////
int shading(Point o, Vector u, Spheres *sph)
{
    while(sph)
    {
          float a = vec_dot(u,u);
          float b = 2.0*(u.x*(o.x - sph->center.x)) + 2.0*(u.y*(o.y - sph->center.y)) + 2.0*(u.z*(o.z - sph->center.z));
          float c = pow(sph->center.x, 2.0) +  pow(sph->center.y, 2.0) + pow(sph->center.z, 2.0) - 
                    2*(sph->center.x*o.x + sph->center.y*o.y + sph->center.z*o.z) + 
                    pow(o.x, 2) +  pow(o.y, 2) + pow(o.z, 2) - 
                    pow(sph->radius, 2);



          float discriminant = pow(b,2) - 4*a*c;
          if(discriminant > 0 )
          {
                  float t0 = (-b - sqrt(discriminant))/2;
                  float t1 = (-b + sqrt(discriminant))/2;
                  if (t0 > 0.0001 || t1 > 0.0001)
                    return 0;// in shadow

          }
          sph = sph->next;
    }
    return 1;

}
/*********************************************************************
 * Phong illumination - you need to implement this!
 *********************************************************************/
 //color = phong(*hit, eye , surf_norm, nearestSph);
RGB_float phong(Point q, Vector v, Vector surf_norm, Spheres *sph, int recursive_step) 
{


    Vector lm = get_vec(q, light1);
    normalize(&lm);

    float distance= vec_len(get_vec(q, light1));

    float angle = vec_dot(surf_norm, vec_scale(lm, -1));
    Vector rm = vec_plus(lm, vec_scale(surf_norm, 2*angle));
    normalize(&rm);


    float decay = 1/(decay_a + decay_b*distance + decay_c*pow(distance,2));

    RGB_float ambient = {0,0,0};
    RGB_float color = {0,0,0};
    ambient.r += global_ambient[0]*sph->reflectance;
    ambient.g += global_ambient[1]*sph->reflectance;
    ambient.b += global_ambient[2]*sph->reflectance;

    ambient.r += global_ambient[0]*sph->reflectance;
    ambient.g += global_ambient[1]*sph->reflectance;
    ambient.b += global_ambient[2]*sph->reflectance;

    color.r += decay*light1_intensity[0]*sph->mat_diffuse[0]*vec_dot(surf_norm, lm);
    color.g += decay*light1_intensity[1]*sph->mat_diffuse[1]*vec_dot(surf_norm, lm); 
    color.b += decay*light1_intensity[2]*sph->mat_diffuse[2]*vec_dot(surf_norm, lm);

    color.r += decay*(light1_intensity[0]*sph->mat_specular[0]*(pow(vec_dot(rm, v), sph->mat_shineness)));
    color.g += decay*(light1_intensity[1]*sph->mat_specular[1]*(pow(vec_dot(rm, v), sph->mat_shineness)));
    color.b += decay*(light1_intensity[2]*sph->mat_specular[2]*(pow(vec_dot(rm, v), sph->mat_shineness)));

    if (shadow_on && (shading(q, lm, scene )==0))
    {
      return ambient;
    }

    if(stochastic_on && recursive_step < max_rays)
    {
       vec3 norm;
       norm.x = surf_norm.x;
       norm.y = surf_norm.y;
       norm.z = surf_norm.z;
       Vector diffuse_dir = vec_scale(v, -1);
       vec3 viewDir;
       viewDir.x = v.x;
       viewDir.y = v.y;
       viewDir.z = v.z;

       vec3 axis = cross(viewDir,norm);
       Vector axis_vector;
       axis_vector.x = axis.x;
       axis_vector.y = axis.y;
       axis_vector.z = axis.z;



       float angle_a = ((float)rand())/10.0 * 20.0 - 10.0;

       float angle_b = ((float)rand())/10.0 * 20.0 - 10.0;

       diffuse_dir = vec_minus(vec_scale(axis_vector, 2*angle_a), diffuse_dir);
       normalize(&diffuse_dir);
       diffuse_dir = vec_minus(vec_scale(surf_norm, 2*angle_b), diffuse_dir);
       normalize(&diffuse_dir);
       RGB_float color_diff = phong(q, diffuse_dir, surf_norm, sph, recursive_step+1);
       color_diff = clr_scale(color_diff, 0.2); 
       color = clr_add(color, color_diff);
    }

    return color;
}


bool isIntersect(Point p, Vector ray, Point *hit)
{
    Vector planeNor = {0,5,3};
    Point planePoint = {-2,-2,-2};

    Vector v;
    v.x = p.x - planePoint.x;
    v.y = p.y - planePoint.y;
    v.z = p.z - planePoint.z;

    float numera = vec_dot(planeNor, v);
    float deno = vec_dot(planeNor, ray);


    if((deno == 0) && (numera != 0)) 
    {
      return false;
    }

    float intersectPoint = -numera/deno;

    if(intersectPoint<=0)
    {
      return false;
    }
    else
    {
      hit->x = p.x + intersectPoint*ray.x;
      hit->y = p.y + intersectPoint*ray.y;
      hit->z = p.z + intersectPoint*ray.z;
      return true;
    }

}

RGB_float getPlaneColor(Point p)
{

  RGB_float black = {0,0,0};
  RGB_float white = {1.0,1.0,1.0};
  if ((int((p.x + 8.0)/2) % 2 != 0 && int((p.z +8.0)/2) % 2 == 0)  || (int((p.x + 8.0)/2) % 2 == 0 && int((p.z +8.0)/2) % 2 != 0) )
  {
      return white;
  }
  return black;

}


/************************************************************************
 * This is the recursive ray tracer - you need to implement this!
 * You should decide what arguments to use.
 ************************************************************************/
RGB_float recursive_ray_trace(Point epos, Vector ray, int step) 
{
    Point *hit = new Point;
    Point *phit = new Point;
    Spheres *nearestSph = intersect_scene(epos, ray, scene, hit, 0);


    RGB_float color = background_clr;


    if(chessboard_on && isIntersect(epos, ray, phit))
    {
        Vector light = get_vec(*phit, epos);
        normalize(&light);
        Vector eye = get_vec(*phit,eye_pos);
        normalize(&eye);
        Vector shadow = get_vec(*phit, light1);

        color = getPlaneColor(*phit);
        if((shadow_on) && (shading(*phit, shadow, scene)==0))
        {
            color = clr_scale(color, 0.5);
        }
       
    }

    if(nearestSph != NULL)
    {
          Vector eye = get_vec(*hit, eye_pos);
          normalize(&eye);
          Vector surf_norm = sphere_normal(*hit, nearestSph);
          normalize(&surf_norm);
          Vector light = get_vec(*hit, epos);
          normalize(&light);

          color = phong(*hit, eye , surf_norm, nearestSph, 0);

          RGB_float reflectionColor;
          RGB_float refractionColor;
          if(reflection_on && step < step_max)
          {

             step++;
             float dot = vec_dot(light, surf_norm);
             Vector scale = vec_scale(surf_norm, 2*dot);
             Vector reflection = vec_minus(scale, light);
             normalize(&reflection);
             reflectionColor = recursive_ray_trace( *hit, reflection, step);
             reflectionColor = clr_scale(reflectionColor, nearestSph->reflectance);
             color = clr_add(color, reflectionColor);

          }

          if(refraction_on && step < step_max)
          {

            Vector view = get_vec(*hit, epos); 
            step++;
            float ratioIndex = nearestSph->refractance/1.0;
            Vector refracRay;
            float cosTheta1 = vec_dot(surf_norm, vec_scale(view, -1));
            float cosTheta2 = sqrt(1.0f-pow(ratioIndex,2)*(1-(pow(cosTheta1,2))));

            if(cosTheta1 > 0.0f)
            {
              refracRay = vec_plus(vec_scale(view, ratioIndex), vec_scale(surf_norm, ratioIndex*cosTheta1-cosTheta2));

            }
            else
            {
              refracRay = vec_minus(vec_scale(view, ratioIndex), vec_scale(surf_norm, ratioIndex*cosTheta1-cosTheta2));
            }

            refracRay = vec_scale(refracRay, -1);
            normalize(&refracRay);
            refracRay.x += hit->x;
            refracRay.y += hit->y;
            refracRay.z += hit->z;
            refractionColor = recursive_ray_trace(*hit, refracRay, step);
            //refractionColor = clr_scale(refractionColor, 0.5);
            color = clr_add(color, refractionColor);
 
          }


          return color;
    }


    return color;
}
      

/*********************************************************************
 * This function traverses all the pixels and cast rays. It calls the
 * recursive ray tracer and assign return color to frame
 *
 * You should not need to change it except for the call to the recursive
 * ray tracer. Feel free to change other parts of the function however,
 * if you must.
 *********************************************************************/
void ray_trace() {
  int i, j;
  float x_grid_size = image_width / float(win_width);
  float y_grid_size = image_height / float(win_height);
  float x_start = -0.5 * image_width;
  float y_start = -0.5 * image_height;
  RGB_float ret_color;
  Point cur_pixel_pos;
  Vector ray;

  // ray is cast through center of pixel
  cur_pixel_pos.x = x_start + 0.5 * x_grid_size;
  cur_pixel_pos.y = y_start + 0.5 * y_grid_size;
  cur_pixel_pos.z = image_plane;

  for (i=0; i<win_height; i++) {
    for (j=0; j<win_width; j++) {
      ray = get_vec(eye_pos, cur_pixel_pos);
      int step = 0;
      ret_color = recursive_ray_trace(eye_pos, ray, step);


      if (super_sampling_on)//add other 4 rays per square pixel
      {

            Vector current_ray;
            RGB_float current_ray_color = {0,0,0};
            Point current_ray_pos = cur_pixel_pos;

            for(float dx = -x_grid_size/4; dx <= x_grid_size; dx+=x_grid_size/2)
            {
              for(float dy = -y_grid_size/4; dy <= y_grid_size; dy+=y_grid_size/2)
              {
                current_ray_pos.x = cur_pixel_pos.x + dx;
                current_ray_pos.y = cur_pixel_pos.y + dy;
                current_ray = get_vec(eye_pos, current_ray_pos);
                ret_color = recursive_ray_trace(eye_pos, current_ray, 0);
                current_ray_color = clr_add(current_ray_color, ret_color); 

              }
            }

            ret_color = clr_scale(current_ray_color, 0.2); 
      }
      //
      // You need to change this!!!
      //
      // ret_color = recursive_ray_trace();
      //ret_color = background_clr; // just background for now

      // Parallel rays can be cast instead using below
      //
      // ray.x = ray.y = 0;
      // ray.z = -1.0;
      // ret_color = recursive_ray_trace(cur_pixel_pos, ray, 1);

      // Checkboard for testing
      // RGB_float clr = {float(i/32), 0, float(j/32)};
      // ret_color = clr;

      frame[i][j][0] = GLfloat(ret_color.r);
      frame[i][j][1] = GLfloat(ret_color.g);
      frame[i][j][2] = GLfloat(ret_color.b);

      cur_pixel_pos.x += x_grid_size;
    }

    cur_pixel_pos.y += y_grid_size;
    cur_pixel_pos.x = x_start;
  }
}
