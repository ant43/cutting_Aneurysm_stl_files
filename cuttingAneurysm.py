from cmath import nan
from re import X
#from tkinter import _XYScrollCommand
import numpy as np
import scipy
import trimesh
from stl import mesh
from mpl_toolkits import mplot3d
from matplotlib import pyplot
from shapely import ops
import rtree
import pyglet
import math
import functools


'''
This section of the code is to detect the Centerline Endpoints
This is one of four steps to creating a cut at the Aneurysm'''


def get_faces(stl_file):
    '''This function will return a numpy array that is all the faces of the .stl file
    'stl_file' this is an stl 3d image file that we will be getting a list of faces from'''
    return trimesh.load_mesh(stl_file, enable_post_processing=True, solid=True).faces


def get_vertices(stl_file):
    '''This function will return a numpy array that is all the vertices of the .stl file
    'stl_file' this is an stl 3d image file that we will be getting a list of vertices from'''
    return trimesh.load_mesh(stl_file, enable_post_processing=True, solid=True).vertices



def get_distance_between_two_points(point1, point2):
    '''point one and point two are two 3 dimentinal points represented by np.arrays of size 3X1 which 
    this function will return the distane between them'''
    return pow(sum(
        [pow(point1[i]-point2[i], 2) for i in range (3)]
    )
        , .5)



def find_closest_vertex(vertices, vertex):
    '''This will look at all the vertexes and get on of the points of intersection then use this point to find the closest vertex 
    between these two points, it will return the location of the vertece that is clocest to this point'''
    #get all the distances
    distances = list(
        map(lambda x : get_distance_between_two_points(x, vertex), vertices))
    #find the minimum distantaces and return the point that produced that
    min_index = 0
    for i in range(len(vertices)):
        if distances[min_index] > distances[i]:
            min_index = i
    return min_index
    #return vertices[min_index]


def find_point_between(point1, point2):
    '''point one and point two are two 3 dimentinal points represented by np.arrays of size 3X1 which 
    this function will return a point that is exactly between them'''
    #this just uses the function of 3d geomotry to caculat the halfway point

    return np.array([((point1[i]+point2[i])/2) 
    for i in range(point1.size)]) , get_distance_between_two_points(point1, point2)






def get_vertices_and_opposit_points(locations, index_ray, ray_origins):
    '''This takes the locations that a ray intersects with and the list of vertexes(ray_origins) and the index_ray and retuerns
    a list of tuples which are the location in the list of vertexes on where the vertex and the oppist point for caculating center
    points is located'''
    #I am using find closest vertex becuase givin a list of vertexes and a vertex and will tell you where the vertex is, 
    # this is why i use
    
    return [(index_ray[i],
    find_closest_vertex(ray_origins, 
    locations[i])) for i in range(index_ray.size)]


def get_center_points(stl_file, parellel_slider = 0.95):
    '''This is to get the points and return a list of them that should be inbetween walls of the vesils, it will first get all the normals
    of vertexes and then creat rays form that and see if the rays intersect, it takes in just the stl file for that, it 
    returns a list of points that are the midpoints of the vessiles, it also returns the radius of the 
    point, how you get the point is output[x][0] and how you get the radius of the point at x is output[x][1]'''


    mesh = trimesh.load_mesh(stl_file, enable_post_processing=True, solid=True)
    ray_origins = mesh.vertices
    ray_directions = mesh.vertex_normals
    ray_directions = ray_directions*(-1)
    locations, index_ray, index_tri = mesh.ray.intersects_location(
        ray_origins=ray_origins,
        ray_directions=ray_directions)

    
    get_vertices_and_opposit_point = get_vertices_and_opposit_points(locations, index_ray, ray_origins)

    get_vertices_and_opposit_point = list(filter
    (lambda x: np.dot(mesh.vertex_normals[x[0]],mesh.vertex_normals[x[1]])>parellel_slider,
    get_vertices_and_opposit_point))


    return list(map(lambda x: find_point_between(ray_origins[x[0]], ray_origins[x[1]]),
     get_vertices_and_opposit_point))

    






    


def calculat_guasian_distrabution(mean, standerd_deviation):
    '''This function takes in the mean and standerd deviation of a guasian 3d distrbution and creates
    it using the equation givin in the reaserch paper'''

    mean = np.transpose([mean])
    left_side = pow(2*math.pi*pow(standerd_deviation, 2), -3/2)
    right_side = np.exp(
        ((-np.transpose([mean])*mean)/(2*pow(standerd_deviation, 2))
    ))
    
    return (left_side*right_side)[0][0][0]



def creat_3d_list(x, y, z):
    '''this creats a 3d list filled with 0's of a specific size'''
    return [[ [0 for col in range(x)] for col in range(y)] for row in range(z)]

def add_guasian_distrabution(accumulation_array, mean, standerd_deviation):
    '''This function takes in mean and standerd_deviation, and takes every point in the accumulation_array, and addes the guasian distrabution 
    into the the accumulation_array, and returns the accumulation array'''
    for x in range(len(accumulation_array)):
        for y in range(len(accumulation_array[0])):
            for z in range(len(accumulation_array[0][0])):
                accumulation_array[x][y][z] = calculat_guasian_distrabution(mean - np.array([x, y, z]), standerd_deviation)
    return accumulation_array

def add_accumulation_arrays(accumulation_array1, accumulation_array2):
    '''this merges two accumulation arrays together and the arrays must be the same size
    in this'''
    for x in range(len(accumulation_array1)):
        for y in range(len(accumulation_array1[0])):
            for z in range(len(accumulation_array1[0][0])):
                accumulation_array1[x][y][z] = accumulation_array1[x][y][z]+accumulation_array2[x][y][z]
    return accumulation_array1



    



def get_accumulator_array(stl_file):
    center_points = get_center_points(stl_file, parellel_slider = 0.95)
    stl_size = trimesh.load_mesh(stl_file, enable_post_processing=True, solid=True).extents
    X = int(stl_size[0])
    Y = int(stl_size[1])
    Z = int(stl_size[2])
    accumulation_arrays = list(map(
        lambda x: add_guasian_distrabution(accumulation_array = creat_3d_list(X, Y, Z) , mean = x[0],
         standerd_deviation = x[1]/2) 
    , center_points))

    return functools.reduce(add_accumulation_arrays, accumulation_arrays)





















#def detect_centerline_endpoints(stl_file):
    '''This function will return a list of 3-dimentional points that are the centerline endpoints of an stl file
    'stl_file' this is an stl 3d image file of the vessiles of which the centerline endpoints will be returned in a list'''
