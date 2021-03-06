"""
Copyright (C) 2016-2017 Andrew Pellegrino

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
from __future__ import division, print_function

import numpy as np
from healpy import ang2pix
from healpy import pixelfunc
from healpy import boundaries
from ctypes import POINTER, c_double, CDLL
import pymetis
import corruscant

class pixel:
    def __init__(self, Nside, id, count):
        self.id = id
        self.count = count
        self.Nside = Nside
        self.neighbors = set(pixelfunc.get_all_neighbours(Nside, id)[::2])
        self.vertecies = set(map(tuple, boundaries(Nside, id).T))

class template:
    def __init__(self, N_fields, theta=None, phi=None, ra=None, dec=None):

        if theta is not None and phi is not None:
            self.theta = np.array(theta)
            self.phi = np.array(phi)

        else:
            self.theta = np.empty(dec.size)
            self.phi = np.empty(ra.size)

            ra = np.require(ra, requirements='AC', dtype='float64')
            dec = np.require(dec, requirements='AC', dtype='float64')

            corruscant.coordlib.radec2sph64(
                    ra.ctypes.data_as(POINTER(c_double)),
                    dec.ctypes.data_as(POINTER(c_double)),
                    self.theta.ctypes.data_as(POINTER(c_double)),
                    self.phi.ctypes.data_as(POINTER(c_double)),
                    ra.size
                                )

        if not (self.phi.shape == self.theta.shape):
            raise ValueError("Theta and phi arrays must be same size")

        self.N_fields = N_fields
        self.Nside = None

        # maps field IDs to color IDs
        self.color_dict = None

        self.filled_pix = None

    def create_fields(self, Nside=1, min_ratio=1.3):
        pixels = ang2pix(Nside,self.theta,self.phi)
        self.Nside = Nside

        # count the data points in each HEALPix pixel
        nbins = 12*Nside*Nside
        bins = np.arange(nbins+1)
        hist, bins = np.histogram(pixels,bins=bins)
        ids = bins[:-1]

        pix_list = [pixel(Nside,id,hist[id]) for id in ids]

        # step to higher resolution if not enough pixels are populated
        self.filled_pix = list(map(lambda p: p.id, filter(lambda p: p.count > 0, pix_list)))
        npix_filled = len(self.filled_pix)
        if npix_filled < self.N_fields:
            print(
                "Nside = {:d} too small for {:d} fields ({:d} pix with data). "
                "Trying Nside = {:d}...".format(Nside, self.N_fields,
                                              npix_filled, 2*Nside)
                 )
            #del pix_list

            return self.create_fields(Nside=2*Nside, min_ratio=min_ratio)
            
        pix_list = sorted(pix_list, key=lambda x: x.id)
        self.pix_list = pix_list

        adj_dict = {}

        xadj_ctr = 0
        xadj = []
        adjncy = []
        eweights = []
        for pix in pix_list:
            xadj.append(xadj_ctr)
            adj_dict[pix.id] = []
            for n in pix.neighbors:
                adj_dict[pix.id].append(n)
                adjncy.append(n)
                xadj_ctr += 1

                if pix.count == 0:
                    eweights.append(1)
                elif pix_list[n].count == 0:
                    eweights.append(1)
                else:
                    eweights.append(1000)

        xadj.append(xadj_ctr)

        self.adj_dict = adj_dict

        vweights = [pix.count for pix in pix_list]

        npix_out, fld_ids = pymetis.part_graph(self.N_fields, xadj=xadj,
                                        adjncy=adjncy, vweights=vweights,
                                        eweights=eweights,
                                        recursive=True)

        # list of field ids parallel to the point list
        fld_ids = np.array(fld_ids)
        pix_ids = fld_ids[np.array(pixels)]
        self.fld_ids = fld_ids

        hist, bins = np.histogram(pix_ids,bins=np.linspace(0,self.N_fields,self.N_fields+1))

        if np.min(hist) == 0:
            print("At least one field has no population. "
                  "Trying Nside = {:d}...".format(2*Nside))
            return self.create_fields(Nside=2*Nside, min_ratio=min_ratio)

        if float(np.min(hist)) / (np.average(hist)) < min_ratio:
            print("The smallest field has less than {:.3f} ".format(min_ratio) +
                  "times the average population. Trying "
                  "Nside = {:d}...".format(2*Nside)
                  )
            return self.create_fields(Nside=2*Nside, min_ratio=min_ratio)

        # TEMP
        return np.array(pix_ids)

    def assign_fields(self, theta=None, phi=None, ra=None, dec=None):
        if theta is not None and phi is not None:
            theta = np.array(theta)
            phi = np.array(phi)

        else:
            theta = np.empty(dec.size)
            phi = np.empty(ra.size)

            ra = np.require(ra, requirements='AC', dtype='float64')
            dec = np.require(dec, requirements='AC', dtype='float64')

            corruscant.coordlib.radec2sph64(
                    ra.ctypes.data_as(POINTER(c_double)),
                    dec.ctypes.data_as(POINTER(c_double)),
                    theta.ctypes.data_as(POINTER(c_double)),
                    phi.ctypes.data_as(POINTER(c_double)),
                    ra.size
                                )

        pix = ang2pix(self.Nside, theta, phi)
        return self.fld_ids[pix]

    def assign_colors(self):
            color_dict = {}
            field_graph = {}

            for key in self.adj_dict:
                f = self.fld_ids[key]
                for entry in self.adj_dict[key]:
                    field_graph.setdefault(self.fld_ids[key], set()).add(self.fld_ids[entry])

            for key in field_graph:
                field_graph[key] = field_graph[key] - set([key])

            # Welsh-Powell algorithm for four-coloring graphs

            color = 0
            sort = sorted(field_graph, key=lambda x: len(field_graph[x]))

            while len(color_dict) < len(field_graph):
                samecolor = []

                for key in reversed(sort):

                    connected = any([neighbor in samecolor for neighbor in field_graph[key]])

                    if not connected:
                        color_dict[key] = color
                        samecolor.append(key)

                for key in color_dict:
                    if key in sort:
                        sort.remove(key)

                color += 1

            return color_dict
        

def main():
    import matplotlib.pyplot as plt

    N_data = 50000

    ######## generate random data ########

    # make it reproducible
    np.random.seed(0)

    # choose uniform points on sphere with independent uniform distributions
    # in cos(theta) and in phi
    costheta_range = np.array([-1, 1])
    phi_range = np.array([0, 2*np.pi])

    costheta = np.random.rand(N_data) * np.diff(costheta_range) + np.min(costheta_range)
    phi = np.random.rand(N_data) * np.diff(phi_range) + np.min(phi_range)

    theta = np.arccos(costheta)

    ######## make fields from data ########

    # maximum acceptable ratio of smallest field count to largest field count
    min_ratio = 0.8

    # number of fields to divide data points into
    N_fields = 20

    templ = template(N_fields, theta, phi)

    # group the dataset into (N_fields) roughly equal groups.
    # Assign a number from 1 to N_fields to each data point
    field_ids = templ.create_fields(min_ratio=min_ratio)
    color_dict = templ.assign_colors()

    hist = np.bincount(field_ids, minlength=N_fields)
    print("from {:d} to {:d}".format(min(field_ids), max(field_ids)))
    print("Least populated field has {:d} points".format(min(hist)))
    print("Most populated field has {:d} points".format(max(hist)))

    ######## plot ########

    fig = plt.figure(figsize=(10,4))
    ax = plt.gca()

    #cm = plt.get_cmap('tab20')

    CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                    '#f781bf', '#a65628', '#984ea3',
                    '#999999', '#e41a1c', '#dede00']
    cm = lambda x: CB_color_cycle[x]

    ax.set_prop_cycle('color',CB_color_cycle)

    field_ids = np.array(field_ids)

    for i in range(N_fields):
        mask = field_ids == i
        phi_field = phi[mask]
        costheta_field = costheta[mask]

        color = cm(color_dict[i])
        plt.plot(phi_field, costheta_field, 'o', markersize=1, color=color)

    plt.title("Jackknife field selection with HEALPix, "
              "n = {:d}".format(N_fields)
              )

    ax.set_xlabel("$\\phi$")
    ax.set_xlim(0, 2 * np.pi)
    ax.set_xticks(np.linspace(0,2*np.pi,5))
    ax.set_xticklabels(["$0$", "$\\pi/2$", "$\\pi$", "$3\\pi/2$", "$2\\pi$"])
    
    ax.set_ylabel("cos $\\theta$")
    ax.set_ylim(-1,1)
    ax.set_yticks(np.linspace(-1,1,5))

    ax.set_aspect('equal')
    plt.tight_layout()
    plt.show()

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
    fig.set_size_inches((10, 4))

    for ax in (ax1, ax2, ax3):
        ax.set_aspect('equal')

    for i in range(N_fields):
        mask = field_ids == i
        color = cm(color_dict[i])

        # negative Y direction
        phi_field = phi[mask & (phi < np.pi)]
        theta_field = theta[mask & (phi < np.pi)]
        ax1.plot(
                np.cos(phi_field)*np.sin(theta_field),
                np.cos(theta_field), 'o', markersize=2,
                color=color
                )

        # positive Y direction
        phi_field = phi[mask & (phi > np.pi)]
        theta_field = theta[mask & (phi > np.pi)]
        ax2.plot(
                np.cos(phi_field)*np.sin(theta_field),
                np.cos(theta_field), 'o', markersize=2,
                color=color
                )

        # negative Z direction
        phi_field = phi[mask & (theta < np.pi/2)]
        theta_field = theta[mask & (theta < np.pi/2)]
        ax3.plot(
                np.cos(phi_field)*np.sin(theta_field),
                np.sin(phi_field)*np.sin(theta_field), 'o', markersize=1,
                color=color
                )


    ax1.set_xlabel("X")
    ax1.invert_xaxis()
    ax1.set_ylabel("Z")

    ax2.set_xlabel("X")
    ax2.set_ylabel("Z")

    ax3.set_xlabel("X")
    ax3.set_ylabel("Y")

    fig.tight_layout()
 
    plt.show()

if __name__ == "__main__":
    main()
