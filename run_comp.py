from sip import init_fields, fields_comp_sip
import json
import os
import analyse_result
import matplotlib.pyplot as plt
import numpy as np
import analyse_result
from db_writer import DB_writer
import time
import datetime
from bson.objectid import ObjectId


def read_params_json(filename):
    if os.path.exists(filename):
        with open(filename, 'r') as json_file:
            return json.load(json_file)

def launch_case(params):
    N = params["N"]
    time_start = datetime.datetime.now()
    psi, T, eta = init_fields(N)
    time_in = 0.
    case_path = analyse_result.generate_case_dir()
    ampl_psi = fields_comp_sip( \
        params["Ra"], \
        params["Le"], \
        params["m"], \
        params["soret"], \
        params["ampl"], \
        params["omega"], \
        time_in, \
        params["time_for_count"], \
        psi, \
        T, \
        eta, \
        case_path
        )
    time_end = datetime.datetime.now()
    delta_time = (time_end - time_start).total_seconds()
    params["ampl_psi"] = ampl_psi
    params["datetime"] = datetime.datetime.utcnow()
    params["comp_time"] = delta_time
    return params

def create_XY_fields(N):
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 1
    x = np.linspace(xmin, xmax, N+1)
    y = np.linspace(ymin, ymax, N+1)
    X, Y = np.meshgrid(x, y)
    return X, Y

if __name__ == "__main__":
    db = DB_writer()
    params = read_params_json('input2.dat')
    for soret in [0., 0.1, 0.2]:
        params["soret"] = soret
        for Le in [0.01, 0.1, 1.0, 10., 100.]:
            params["Le"] = Le
            for Ra in [90, 100, 110, 120]:
                params["Ra"] = Ra
                params = launch_case(params)
                params['_id'] = ObjectId() 
                db.insert_document(params)
                print(soret, Le, Ra)
    db.show_content()
        # TODO: read from file to picture
        # analyse_result.savefig_case()
        # fig, axes = plt.subplots(nrows = 1, ncols = 2)
        # X, Y = create_XY_fields(params["N"])
        # analyse_result.create_field_plot(X, Y, psi.transpose(), u'$\psi$', axes.flat[0])
        # analyse_result.create_field_plot(X, Y, T.transpose(), u'$T$', axes.flat[1])
        # fig.savefig("pic.png")