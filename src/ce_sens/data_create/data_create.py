from ce_sens.data_create.sample import create_data
import argparse

def data_create():
    parser = argparse.ArgumentParser()
    parser.add_argument("inp_path", type=str, help="path to ini file")
    parser.add_argument("rho", help="local merger rate")
    parser.add_argument("time", type=float, help="observation time")
    parser.add_argument("z_max", type=float, help="maximum z")
    parser.add_argument("out_path", help="Output path")

    args = parser.parse_args()
    inp_path = [args.inp_path]
    out_path = args.out_path
    rho = eval(args.rho)
    time = args.time
    z_max = args.z_max

    create_data(inp_path, out_path, rho, time, z_max)
