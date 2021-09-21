import os, shutil, time
import subprocess

EXECUTABLE = 'g09.exe'
GAUSS_EXEDIR = 'C:\\Program Files\\G09W'
SHARED_PROCS = 4

class CalculationError(Exception):
    pass

class WriteParametersError(Exception):
    pass

def start_calc(
    filename,
    input_dir = None,
    output_dir = None,
    scr_dir = None,
    working_dir = None,
    override = False,
    clean = False
):
    calc_name, _ = os.path.splitext(os.path.basename(filename))

    if not input_dir:
        input_dir = os.path.dirname(filename)
    if not output_dir:
        output_dir = os.path.dirname(filename)
    input_dir = os.path.abspath(input_dir)
    output_dir = os.path.abspath(output_dir)

    if working_dir:
        working_dir = os.path.abspath(working_dir)

    runfile = os.path.join(input_dir, os.path.basename(filename))
    outfile = os.path.join(output_dir, calc_name + '.out')

    if not os.path.isfile(runfile):
        raise CalculationError('No file %s' % runfile, runfile)
    if os.path.isfile(outfile) and not override:
        raise CalculationError('Output file %s already exists' % outfile, outfile)

    run_env = os.environ.copy()
    run_env['GAUSS_EXEDIR'] = GAUSS_EXEDIR
    if scr_dir:
        run_env['GAUSS_SCRDIR'] = scr_dir

    return_code = subprocess.run([
        EXECUTABLE,
        runfile,
        outfile
    ], cwd=working_dir, env=run_env).returncode
    if return_code:
        raise CalculationError('Gaussian calculation finished with code %s' % return_code, return_code)
    if working_dir and clean:
        clean_directory(working_dir)

def clean_directory(dir):
    for filename in os.listdir(dir):
        file = os.path.join(dir, filename)
        try:
            if os.path.isfile(file):
                os.unlink(file)
            elif os.path.isdir(file):
                shutil.rmtree(file)
        except OSError as e:
            print('Failed to delete file %s because of %s' % (file, e))

def multiple_calc(
    filenames,
    input_dir = None,
    output_dir = None,
    scr_dir = None,
    working_dir = None,
    override = False,
    clean = False
):
    if input_dir and (not os.path.isdir(input_dir)):
        raise CalculationError('Input folder %s does not exist' % input_dir, input_dir)

    if output_dir and (not os.path.isdir(output_dir)):
        os.makedirs(output_dir)
    if scr_dir and (not os.path.isdir(scr_dir)):
        os.makedirs(scr_dir)
    if working_dir and (not os.path.isdir(working_dir)):
        os.makedirs(working_dir)

    for filename in filenames:
        start_calc(
            filename = filename, 
            input_dir = input_dir,
            output_dir = output_dir,
            scr_dir = scr_dir,
            working_dir = working_dir,
            override = override,
            clean = clean)

def write_input_params(
    coord_filename,
    out_filename = None,
    method = 'b3lyp',
    basis = '6-311++g',
    solvent = None,
    job = 'opt',
    nprocs = SHARED_PROCS,
    coord_dir = None,
    output_dir = None,
    chk_dir = None
):
    calc_name, _ = os.path.splitext(os.path.basename(out_filename or coord_filename))

    if coord_dir:
        coordfile = os.path.join(coord_dir, coord_filename)
    else:
        coordfile = coord_filename
    coordfile = os.path.abspath(coordfile)

    if output_dir and out_filename:
        outfile = os.path.join(output_dir, out_filename)
    elif output_dir:
        outfile = os.path.join(output_dir, calc_name + '.gjf')
    elif out_filename:
        outfile = out_filename
    else:
        outfile = calc_name + '.gjf'
    outfile = os.path.abspath(outfile)

    if chk_dir:
        chkfile = os.path.abspath(os.path.join(chk_dir, calc_name + '.chk'))

    coords = None
    if coordfile == outfile:
        cf = open(coordfile)
        for line in cf:
            coords += line
        cf.close()

    if not os.path.isfile(coordfile):
        raise WriteParametersError('No file %s' % coordfile, coordfile)
    if solvent:
        solvent = 'scrf=(cpcm,solvent=%s)' % solvent

    of = open(outfile, 'w+')
    if chk_dir:
        of.write('%%chk=%s\n' % chkfile)
    of.write('%%nprocshared=%i\n' % nprocs)
    of.write('# %s %s %s %s\n\n' % (
        method,
        basis,
        solvent,
        job
    ))
    of.write('Title: %s %s\n\n' % (calc_name, job))

    if coords:
        of.write(coords)
    else:
        cf = open(coordfile)
        for line in cf:
            of.write(line)
        cf.close()

    of.write('\n\n\n')
    of.close()

def write_multiple_input_params(
    coord_filenames,
    method = 'b3lyp',
    basis = '6-311++g',
    solvent = None,
    job = 'opt',
    nprocs = SHARED_PROCS,
    coord_dir = None,
    output_dir = None,
    chk_dir = None
):
    if coord_dir and (not os.path.isdir(coord_dir)):
        raise WriteParametersError('Coordinates folder %s does not exist' % coord_dir, coord_dir)
    
    if output_dir and (not os.path.isdir(output_dir)):
        os.makedirs(output_dir)

    for coord_filename in coord_filenames:
        write_input_params(
            coord_filename = coord_filename,
            method = method,
            basis = basis,
            solvent = solvent,
            job = job,
            nprocs = nprocs,
            coord_dir = coord_dir,
            output_dir = output_dir,
            chk_dir = chk_dir
        )