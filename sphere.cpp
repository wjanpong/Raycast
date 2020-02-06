#include "sphere.h"
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <cstdlib>
#include <iostream>
/**********************************************************************
 * This function intersects a ray with a given sphere 'sph'. You should
 * use the parametric representation of a line and do the intersection.
 * The function should return the parameter value for the intersection, 
 * which will be compared with others to determine which intersection
 * is closest. The value -1.0 is returned if there is no intersection
 *
 * If there is an intersection, the point of intersection should be
 * stored in the "hit" variable
 **********************************************************************/
float intersect_sphere(Point o, Vector u, Spheres *sph, Point *hit) 
{
	  float a = vec_dot(u,u);
    float b = 2.0*(u.x*(o.x - sph->center.x)) + 2.0*(u.y*(o.y - sph->center.y)) + 2.0*(u.z*(o.z - sph->center.z));
    float c = pow(sph->center.x, 2.0) +  pow(sph->center.y, 2.0) + pow(sph->center.z, 2.0) - 
              2*(sph->center.x*o.x + sph->center.y*o.y + sph->center.z*o.z) + 
              pow(o.x, 2) +  pow(o.y, 2) + pow(o.z, 2) - 
              pow(sph->radius, 2);



    float discriminant = pow(b,2) - 4*a*c;
    if(discriminant<0)
    {
      return -1.0;
    }
    else
    {
          float t = (-b - sqrt(discriminant))/(2.0*a);

          if(t<0)
            return -1.0;

          hit->x = o.x + t*u.x;
          hit->y = o.y + t*u.y;
          hit->z = o.z + t*u.z;

          Vector intersect =  { hit->x-o.x , hit->y-o.y , hit->z-o.z};
          float intersectDis = vec_len(intersect);
          return intersectDis;

    }
}



/*********************************************************************
 * This function returns a pointer to the sphere object that the
 * ray intersects first; NULL if no intersection. You should decide
 * which arguments to use for the function. For exmaple, note that you
 * should return the point of intersection to the calling function.
 **********************************************************************/
Spheres *intersect_scene(Point o, Vector u, Spheres *sphList, Point *hit, int value) 
{
//
// do your thing here
//
    float nearestIntersect = FLT_MAX;
    float intersectDistant;
    Spheres *nearestSph = NULL;
  	Spheres *currentSph = sphList;

    while (currentSph != NULL)
    {
      intersectDistant = intersect_sphere(o, u, currentSph, hit);
      if((intersectDistant>=0) && (intersectDistant < nearestIntersect))
      {
        nearestIntersect = intersectDistant;
        nearestSph = currentSph;
      }
      currentSph = currentSph->next;
    }

    return nearestSph;
}

Spheres *add_sphere(Spheres *slist, Point ctr, float rad, float amb[],
        float dif[], float spe[], float shine, 
        float refl, int sindex, float refra) {
  Spheres *new_sphere;

  new_sphere = (Spheres *)malloc(sizeof(Spheres));
  new_sphere->index = sindex;
  new_sphere->center = ctr;
  new_sphere->radius = rad;
  (new_sphere->mat_ambient)[0] = amb[0];
  (new_sphere->mat_ambient)[1] = amb[1];
  (new_sphere->mat_ambient)[2] = amb[2];
  (new_sphere->mat_diffuse)[0] = dif[0];
  (new_sphere->mat_diffuse)[1] = dif[1];
  (new_sphere->mat_diffuse)[2] = dif[2];
  (new_sphere->mat_specular)[0] = spe[0];
  (new_sphere->mat_specular)[1] = spe[1];
  (new_sphere->mat_specular)[2] = spe[2];
  new_sphere->mat_shineness = shine;
  new_sphere->reflectance = refl;
  new_sphere->refractance = refra;
  new_sphere->next = NULL;

  if (slist == NULL) { // first object
    slist = new_sphere;
  } else { // insert at the beginning
    new_sphere->next = slist;
    slist = new_sphere;
  }

  return slist;
}


/*****************************************************
 * This function adds a sphere into the sphere list
 *
 * You need not change this.
 *****************************************************/
Spheres *add_sphere(Spheres *slist, Point ctr, float rad, float amb[],
		    float dif[], float spe[], float shine, 
		    float refl, int sindex) {
  Spheres *new_sphere;

  new_sphere = (Spheres *)malloc(sizeof(Spheres));
  new_sphere->index = sindex;
  new_sphere->center = ctr;
  new_sphere->radius = rad;
  (new_sphere->mat_ambient)[0] = amb[0];
  (new_sphere->mat_ambient)[1] = amb[1];
  (new_sphere->mat_ambient)[2] = amb[2];
  (new_sphere->mat_diffuse)[0] = dif[0];
  (new_sphere->mat_diffuse)[1] = dif[1];
  (new_sphere->mat_diffuse)[2] = dif[2];
  (new_sphere->mat_specular)[0] = spe[0];
  (new_sphere->mat_specular)[1] = spe[1];
  (new_sphere->mat_specular)[2] = spe[2];
  new_sphere->mat_shineness = shine;
  new_sphere->reflectance = refl;
  new_sphere->next = NULL;

  if (slist == NULL) { // first object
    slist = new_sphere;
  } else { // insert at the beginning
    new_sphere->next = slist;
    slist = new_sphere;
  }

  return slist;
}

/******************************************
 * computes a sphere normal - done for you
 ******************************************/
Vector sphere_normal(Point q, Spheres *sph) {
  Vector rc;

  rc = get_vec(sph->center, q);
  normalize(&rc);
  return rc;
}
