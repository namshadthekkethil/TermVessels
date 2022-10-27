// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

// <h1> Systems Example 7 - Large deformation elasticity (St. Venant-Kirchoff
// material) </h1> \author Lorenzo Zanon \author David Knezevic \date 2014
//
// In this example, we consider an elastic cantilever beam modeled as a St.
// Venant-Kirchoff material (which is an extension of the linear elastic
// material model to the nonlinear regime). The implementation presented here
// uses NonlinearImplicitSystem.
//
// We formulate the PDE on the reference geometry (\Omega) as opposed to the
// deformed geometry (\Omega^deformed). As a result (e.g. see Ciarlet's 3D
// elasticity book, Theorem 2.6-2) the PDE is given as follows:
//
//     \int_\Omega F_im Sigma_mj v_i,j = \int_\Omega f_i v_i + \int_\Gamma g_i
//     v_i ds
//
// where:
//  * F is the deformation gradient, F = I + du/dx (x here refers to reference
//  coordinates).
//  * Sigma is the second Piola-Kirchoff stress, which for the St. Venant
//  Kirchoff model is
//    given by Sigma_ij = C_ijkl E_kl, where E_kl is the strain,
//    E_kl = 0.5 * (u_k,l + u_l,k + u_m,k u_m,l).
//  * f is a body load.
//  * g is a surface traction on the surface \Gamma.
//
// In this example we only consider a body load (e.g. gravity), hence we set g =
// 0.

// C++ include files that we need
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>

// Various include files needed for the mesh & solver functionality.
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/equation_systems.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/fe.h"
#include "libmesh/getpot.h"
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/zero_function.h"

#include "libmesh/boundary_info.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/edge_edge3.h"
#include "libmesh/point.h"

#include "libmesh/petsc_matrix.h"

#include "libmesh/tecplot_io.h"

#include "libmesh/linear_implicit_system.h"

// The nonlinear solver and system we will be using
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/nonlinear_solver.h"

#include <libmesh/mesh_modification.h>
#include <libmesh/mesh_refinement.h>
#include <libmesh/mesh_tools.h>

#include "libmesh/mesh_tetgen_interface.h"

using namespace libMesh;
using namespace std::chrono;

using namespace std;

typedef struct Vess {
  int p1, p2, p, dl, dr, i_g;
  double r;
} Vess;

vector<Vess> vessels;
Vess vess_i;
int i_group, i_g_count;
vector<double> r_group;

void inside_elements(Mesh &mesh, double xl, double yl, double zl,
                     int &inside_LV) {
  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
      mesh.active_local_elements_end();

  inside_LV = 0;

  for (; el != end_el; ++el) {
    Elem *elem = *el;
    const int elem_id = elem->id();

    double x_elem = 0.0;
    double y_elem = 0.0;
    double z_elem = 0.0;
    for (unsigned int j = 0; j < 4; j++) {

      Point pj;
      pj = elem->point(j);

      x_elem += pj(0);
      y_elem += pj(1);
      z_elem += pj(2);
    }
    x_elem /= 4.0;
    y_elem /= 4.0;
    z_elem /= 4.0;

    double h_elem = min(elem->length(0, 1), elem->length(1, 2));
    h_elem = min(h_elem, elem->length(0, 2));

    h_elem = min(h_elem, elem->length(0, 3));
    h_elem = min(h_elem, elem->length(1, 3));
    h_elem = min(h_elem, elem->length(2, 3));

    double dist_elem =
        sqrt(pow(xl - x_elem, 2) + pow(yl - y_elem, 2) + pow(zl - z_elem, 2));

    if (dist_elem < 0.5 * h_elem)
      inside_LV = 1;
  }
}

void update_vessel(vector<int> &p1, vector<int> &p2, vector<double> &x,
                   vector<double> &y, vector<double> &z, vector<int> &dl,
                   vector<int> &dr, vector<int> &pr, int i, Mesh &mesh) {
  // cout<<"i="<<i<<" dl[i]="<<dl[i]<<endl;
  if (dl[i] != -10) {
    if (dr[i] == -10) {
      vess_i.p1 = p1[dl[i]];
      vess_i.p2 = p2[dl[i]];
      vess_i.p = -10;
      vess_i.dl = -10;
      vess_i.dr = -10;

      vessels.push_back(vess_i);
    }

    else {
      vess_i.p1 = p1[dl[i]];
      vess_i.p2 = p2[dl[i]];
      vess_i.p = -10;
      vess_i.dl = -10;
      vess_i.dr = -10;

      vessels.push_back(vess_i);

      vess_i.p1 = p1[dr[i]];
      vess_i.p2 = p2[dr[i]];
      vess_i.p = -10;
      vess_i.dl = -10;
      vess_i.dr = -10;

      vessels.push_back(vess_i);
    }

    double xl, yl, zl;
    xl = 0.5 * (x[p1[i]] + x[p2[i]]);
    yl = 0.5 * (y[p1[i]] + y[p2[i]]);
    zl = 0.5 * (z[p1[i]] + z[p2[i]]);

    int inside_LV = 0;

    // if (dr[pr[i]] != -10)
    //   inside_elements(mesh, xl, yl, zl, inside_LV);
    //
    // if (inside_LV == 0) {
    //   update_vessel(p1, p2, x, y, z, dl, dr, pr, dl[i], mesh);
    //
    //   if (dr[i] != -10)
    //     update_vessel(p1, p2, x, y, z, dl, dr, pr, dr[i], mesh);
    // }

    if (dr[dl[i]] != -10)
      inside_elements(mesh, xl, yl, zl, inside_LV);

    if (inside_LV == 0)
      update_vessel(p1, p2, x, y, z, dl, dr, pr, dl[i], mesh);

    if (dr[i] != -10) {
      inside_LV = 0;
      if (dr[dr[i]] != -10)
        inside_elements(mesh, xl, yl, zl, inside_LV);

      if (inside_LV == 0)
        update_vessel(p1, p2, x, y, z, dl, dr, pr, dr[i], mesh);
    }
  }
}

void update_parent_daughter() {
  for (int i = 0; i < vessels.size(); i++) {
    vessels[i].p = -10;
    int count_pr = 0;
    for (int j = 0; j < vessels.size(); j++) {
      if (vessels[i].p1 == vessels[j].p2) {
        vessels[i].p = j;
        count_pr++;
      }
    }

    // if(count_pr > 1)
    // pr[i] = 1000000;

    if (count_pr > 1 || count_pr < 1)
      cout << "i=" << i << " count_pr=" << count_pr
           << " error assigning parents" << endl;
  }

  for (int i = 0; i < vessels.size(); i++) {
    if (vessels[i].p != -10) {
      if (vessels[vessels[i].p].dl == -10)
        vessels[vessels[i].p].dl = i;
      else if (vessels[vessels[i].p].dr == -10)
        vessels[vessels[i].p].dr = i;
      else
        cout << "i=" << i << " error assigning daughters" << endl;
    }

    // if(pr[i]== -10)
    // cout<<"i="<<i<<" pr="<<pr[i]<<endl;
    // else
    // cout<<"i="<<i<<" pr="<<pr[i]<<" dl="<<dl[pr[i]]<<" dr="<<dr[pr[i]]<<endl;
  }
}

void update_radius(int n) {
  if (vessels[n].dl != -10 && vessels[n].dr == -10) {
    r_group[i_group] += vessels[vessels[n].dl].r;
    vessels[vessels[n].dl].i_g = i_group;
    i_g_count++;

    update_radius(vessels[n].dl);
  } else {
    r_group[i_group] /= i_g_count;
    i_group++;
    i_g_count = 0;

    r_group.push_back(0.0);

    if (vessels[n].dr != -10) {
      r_group[i_group] += vessels[vessels[n].dl].r;
      vessels[vessels[n].dl].i_g = i_group;
      i_g_count++;

      update_radius(vessels[n].dl);

      r_group[i_group] /= i_g_count;
      i_group++;
      i_g_count = 0;
      r_group.push_back(0.0);

      r_group[i_group] += vessels[vessels[n].dr].r;
      vessels[vessels[n].dr].i_g = i_group;
      i_g_count++;

      update_radius(vessels[n].dr);
    }
  }
}

void assign_radius(vector<double> &r, vector<int> &p1, vector<int> &p2) {
  for (int i = 0; i < vessels.size(); i++) {
    vessels[i].r = 0.5 * (r[vessels[i].p1] + r[vessels[i].p2]);
  }

  i_group = 0;
  r_group.push_back(0.0);

  i_g_count = 0;

  r_group[i_group] += vessels[0].r;
  vessels[0].i_g = i_group;
  i_g_count++;
  update_radius(0);

  for (int i = 0; i < vessels.size(); i++) {
    vessels[i].r = r_group[vessels[i].i_g];
  }
}

int main(int argc, char **argv) {

  LibMeshInit init(argc, argv);

  // This example requires the PETSc nonlinear solvers
  // libmesh_example_requires(libMesh::default_solver_package() ==
  // PETSC_SOLVERS,
  //                          "--enable-petsc");

  // We use a 3D domain.
  libmesh_example_requires(LIBMESH_DIM > 2,
                           "--disable-1D-only --disable-2D-only");

  int rank, np;

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  Mesh mesh(init.comm());

  mesh.read("LVtrunk_heart_real_Lee_coupl.xda", NULL);

  vector<double> x, y, z, r;
  vector<int> p1, p2, pr, dl, dr, d3, d4, p1_new, p2_new;
  vector<double> x1, y1, z1, x2, y2, z2, l, r1, r2;

  p1.resize(17903);
  p2.resize(17903);
  ifstream file_p;
  file_p.open("edges.dat");
  for (int i = 0; i < 17903; i++) {
    file_p >> p1[i] >> p2[i];
  }
  file_p.close();

  for (int i = 0; i < 17903; i++) {
    p1[i] -= 1;
    p2[i] -= 1;
  }

  x.resize(17904);
  y.resize(17904);
  z.resize(17904);
  ifstream file_x;
  file_x.open("vasc.dat");
  for (int i = 0; i < 17904; i++) {
    file_x >> x[i] >> y[i] >> z[i];
  }
  file_x.close();

  r.resize(17904);
  ifstream file_r;
  file_r.open("rad.dat");
  for (int i = 0; i < 17904; i++) {
    file_r >> r[i];
    // r[i] = 1.0;//
  }
  file_r.close();

  // for(int i=0;i<17903;i++)
  //{
  // x1.push_back(x[p1[i]]);
  // y1.push_back(y[p1[i]]);
  // z1.push_back(z[p1[i]]);

  // x2.push_back(x[p2[i]]);
  // y2.push_back(y[p2[i]]);
  // z2.push_back(z[p2[i]]);

  // l.push_back(sqrt((x2[i]-x1[i])*(x2[i]-x1[i])+(y2[i]-y1[i])*(y2[i]-y1[i])+(z2[i]-z1[i])*(z2[i]-z1[i])));

  // r1.push_back(r[p1[i]]);
  // r2.push_back(r[p2[i]]);

  //}

  // int i=17903;

  p1.push_back(16703);
  p2.push_back(16714);

  // x1.push_back(x[p1[i]]);
  // y1.push_back(y[p1[i]]);
  // z1.push_back(z[p1[i]]);

  // x2.push_back(x[p2[i]]);
  // y2.push_back(y[p2[i]]);
  // z2.push_back(z[p2[i]]);

  // l.push_back(sqrt((x2[i]-x1[i])*(x2[i]-x1[i])+(y2[i]-y1[i])*(y2[i]-y1[i])+(z2[i]-z1[i])*(z2[i]-z1[i])));

  // r1.push_back(r[p1[i]]);
  // r2.push_back(r[p2[i]]);

  // i=17904;

  p1.push_back(8852);
  p2.push_back(8926);

  // x1.push_back(x[p1[i]]);
  // y1.push_back(y[p1[i]]);
  // z1.push_back(z[p1[i]]);

  // x2.push_back(x[p2[i]]);
  // y2.push_back(y[p2[i]]);
  // z2.push_back(z[p2[i]]);

  // l.push_back(sqrt((x2[i]-x1[i])*(x2[i]-x1[i])+(y2[i]-y1[i])*(y2[i]-y1[i])+(z2[i]-z1[i])*(z2[i]-z1[i])));

  // r1.push_back(r[p1[i]]);
  // r2.push_back(r[p2[i]]);

  int p1_temp = p2[15158];
  int p2_temp = p1[15158];

  p1[15158] = p1_temp;
  p2[15158] = p2_temp;

  p1_temp = p2[14843];
  p2_temp = p1[14843];

  p1[14843] = p1_temp;
  p2[14843] = p2_temp;

  p1_temp = p2[9675];
  p2_temp = p1[9675];

  p1[9675] = p1_temp;
  p2[9675] = p2_temp;

  int i_new = 17904;
  for (int ii = 0; ii < 17905; ii++) {
    // double x_mid = 0.5*(x[p1[ii]]+x[p2[ii]]);
    // double y_mid = 0.5*(y[p1[ii]]+y[p2[ii]]);
    // double z_mid = 0.5*(z[p1[ii]]+z[p2[ii]]);
    // double r_mid = 0.5*(r[p1[ii]]+r[p2[ii]]);

    // x.push_back(x_mid);
    // y.push_back(y_mid);
    // z.push_back(z_mid);
    // r.push_back(r_mid);

    // p1_new.push_back(p1[ii]);
    // p2_new.push_back(i_new);

    // p1_new.push_back(i_new);
    // p2_new.push_back(p2[ii]);

    // i_new++;

    double x_1 = 0.75 * x[p1[ii]] + 0.25 * x[p2[ii]];
    double y_1 = 0.75 * y[p1[ii]] + 0.25 * y[p2[ii]];
    double z_1 = 0.75 * z[p1[ii]] + 0.25 * z[p2[ii]];
    double r_1 = 0.75 * r[p1[ii]] + 0.25 * r[p2[ii]];

    double x_2 = 0.5 * (x[p1[ii]] + x[p2[ii]]);
    double y_2 = 0.5 * (y[p1[ii]] + y[p2[ii]]);
    double z_2 = 0.5 * (z[p1[ii]] + z[p2[ii]]);
    double r_2 = 0.5 * (r[p1[ii]] + r[p2[ii]]);

    double x_3 = 0.25 * x[p1[ii]] + 0.75 * x[p2[ii]];
    double y_3 = 0.25 * y[p1[ii]] + 0.75 * y[p2[ii]];
    double z_3 = 0.25 * z[p1[ii]] + 0.75 * z[p2[ii]];
    double r_3 = 0.25 * r[p1[ii]] + 0.75 * r[p2[ii]];

    x.push_back(x_1);
    y.push_back(y_1);
    z.push_back(z_1);
    r.push_back(r_1);

    p1_new.push_back(p1[ii]);
    p2_new.push_back(i_new);

    i_new++;

    x.push_back(x_2);
    y.push_back(y_2);
    z.push_back(z_2);
    r.push_back(r_2);

    p1_new.push_back(i_new - 1);
    p2_new.push_back(i_new);

    i_new++;

    x.push_back(x_3);
    y.push_back(y_3);
    z.push_back(z_3);
    r.push_back(r_3);

    p1_new.push_back(i_new - 1);
    p2_new.push_back(i_new);

    p1_new.push_back(i_new);
    p2_new.push_back(p2[ii]);

    i_new++;
  }

  p1.resize(p1_new.size());
  p2.resize(p1_new.size());

  for (int i = 0; i < p1_new.size(); i++) {
    p1[i] = p1_new[i];
    p2[i] = p2_new[i];
  }

  for (int i = 0; i < p1.size(); i++) {
    x1.push_back(x[p1[i]]);
    y1.push_back(y[p1[i]]);
    z1.push_back(z[p1[i]]);

    x2.push_back(x[p2[i]]);
    y2.push_back(y[p2[i]]);
    z2.push_back(z[p2[i]]);

    l.push_back(sqrt((x2[i] - x1[i]) * (x2[i] - x1[i]) +
                     (y2[i] - y1[i]) * (y2[i] - y1[i]) +
                     (z2[i] - z1[i]) * (z2[i] - z1[i])));

    r1.push_back(r[p1[i]]);
    r2.push_back(r[p2[i]]);
  }

  pr.resize(l.size());
  for (int i = 0; i < l.size(); i++) {
    pr[i] = -10;
    int count_pr = 0;
    for (int j = 0; j < l.size(); j++) {
      if (p1[i] == p2[j]) {
        pr[i] = j;
        count_pr++;
      }
    }

    // if(count_pr > 1)
    // pr[i] = 1000000;

    if (count_pr > 1 || count_pr < 1)
      cout << "i=" << i << " count_pr=" << count_pr
           << " error assigning parents" << endl;
  }

  dl.resize(l.size());
  dr.resize(l.size());
  d3.resize(l.size());
  d4.resize(l.size());
  for (int i = 0; i < l.size(); i++) {
    dl[i] = -10;
    dr[i] = -10;
    d3[i] = -10;
    d4[i] = -10;
  }
  for (int i = 0; i < l.size(); i++) {
    if (pr[i] != -10) {
      if (dl[pr[i]] == -10)
        dl[pr[i]] = i;
      else if (dr[pr[i]] == -10)
        dr[pr[i]] = i;
      else if (d3[pr[i]] == -10)
        d3[pr[i]] = i;
      else if (d4[pr[i]] == -10)
        d4[pr[i]] = i;
      else
        cout << "i=" << i << " error assigning daughters" << endl;
    }

    // if(pr[i]== -10)
    // cout<<"i="<<i<<" pr="<<pr[i]<<endl;
    // else
    // cout<<"i="<<i<<" pr="<<pr[i]<<" dl="<<dl[pr[i]]<<" dr="<<dr[pr[i]]<<endl;
  }

  ofstream file1;
  file1.open("terminals_Lee_full.dat", ios::out);
  for (int i = 0; i < l.size(); i++) {
    int count_pri = 0;
    if (pr[i] != -10) {
      for (int j = 0; j < l.size(); j++) {
        if (pr[i] == pr[j])
          count_pri++;
      }
    }

    file1 << x1[i] << "," << y1[i] << "," << z1[i] << "," << x2[i] << ","
          << y2[i] << "," << z2[i] << "," << l[i] << "," << r1[i] << ","
          << r2[i] << "," << pr[i] << "," << p1[i] << "," << p2[i] << "," << i
          << "," << count_pri << "," << dl[i] << "," << dr[i] << "," << d3[i]
          << "," << d4[i] << endl;
  }
  file1.close();

  ofstream file3;
  file3.open("terminals_Lee_ter.dat", ios::out);
  for (int i = 0; i < l.size(); i++) {
    if (dl[i] == -10)
      file3 << x1[i] << " " << y1[i] << " " << z1[i] << " " << x2[i] << " "
            << y2[i] << " " << z2[i] << " " << l[i] << " "
            << 0.5 * (r1[i] + r2[i]) << endl;
  }
  file3.close();

  int count_st_vess = 0;
  int start_vess;
  for (int i = 0; i < l.size(); i++) {
    if (pr[i] == -10) {
      start_vess = i;
      count_st_vess++;
    }
  }
  // if(count_st_vess != 1)
  // cout<<"i="<<i<<" error assigning daughters"<<endl;

  vess_i.p1 = p1[start_vess];
  vess_i.p2 = p2[start_vess];
  vess_i.p = -10;
  vess_i.dl = -10;
  vess_i.dr = -10;

  vessels.push_back(vess_i);

  update_vessel(p1, p2, x, y, z, dl, dr, pr, start_vess, mesh);

  update_parent_daughter();

  assign_radius(r, p1, p2);

  ofstream file_vess;
  file_vess.open("updated_vessels_Lee_term.csv", ios::out);
  file_vess << "\"x1\""
            << ",\"y1\""
            << ",\"z1\""
            << ",\"x2\""
            << ",\"y2\""
            << ",\"z2\""
            << ",\"l\""
            << ",\"r1\""
            << ",\"r2\""
            << ",\"pr\""
            << ",\"dl\""
            << ",\"dr\""
            << ",\"inside\"" << endl;
  for (int i = 0; i < vessels.size(); i++) {
    double x1_vess = x[vessels[i].p1];
    double y1_vess = y[vessels[i].p1];
    double z1_vess = z[vessels[i].p1];

    double x2_vess = x[vessels[i].p2];
    double y2_vess = y[vessels[i].p2];
    double z2_vess = z[vessels[i].p2];

    double l_vess = sqrt((x2_vess - x1_vess) * (x2_vess - x1_vess) +
                         (y2_vess - y1_vess) * (y2_vess - y1_vess) +
                         (z2_vess - z1_vess) * (z2_vess - z1_vess));
    double r1_vess = vessels[i].r; // r[vessels[i].p1];
    double r2_vess = vessels[i].r; // r[vessels[i].p2];

    file_vess << x1_vess << "," << y1_vess << "," << z1_vess << "," << x2_vess
              << "," << y2_vess << "," << z2_vess << "," << l_vess << ","
              << r1_vess << "," << r2_vess << "," << vessels[i].p << ","
              << vessels[i].dl << "," << vessels[i].dr << "," << 0 << endl;

    // cout<<"i="<<i<<" r1="<<r1_vess<<" r2="<<r2_vess<<"
    // "<<r1_vess-r2_vess<<endl;
  }
  file_vess.close();

  ofstream file_vess_dat;
  file_vess_dat.open("updated_vessels_Lee_data_term.dat", ios::out);

  for (int i = 0; i < vessels.size(); i++) {
    double x1_vess = x[vessels[i].p1];
    double y1_vess = y[vessels[i].p1];
    double z1_vess = z[vessels[i].p1];

    double x2_vess = x[vessels[i].p2];
    double y2_vess = y[vessels[i].p2];
    double z2_vess = z[vessels[i].p2];

    double l_vess = sqrt((x2_vess - x1_vess) * (x2_vess - x1_vess) +
                         (y2_vess - y1_vess) * (y2_vess - y1_vess) +
                         (z2_vess - z1_vess) * (z2_vess - z1_vess));
    double r1_vess = vessels[i].r; // r[vessels[i].p1];
    double r2_vess = vessels[i].r; // r[vessels[i].p2];
    file_vess_dat << x1_vess << " " << y1_vess << " " << z1_vess << " "
                  << x2_vess << " " << y2_vess << " " << z2_vess << " "
                  << l_vess << " " << r1_vess << " " << r2_vess << " "
                  << 0.5 * (r1_vess + r2_vess) << " " << vessels[i].p << " "
                  << vessels[i].dl << " " << vessels[i].dr << " " << 0 << endl;

    // cout<<"i="<<i<<" r1="<<r1_vess<<" r2="<<r2_vess<<"
    // "<<r1_vess-r2_vess<<endl;
  }
  file_vess_dat.close();
}
