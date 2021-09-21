from shutil import ExecError
import numpy as np

LINE_SEPARATOR = '-' * 20
PARTS_SEPARATOR = '*' * 20
HARTREE_TO_KCAL_PER_MOLE = 627.5

ELEMENTS = {
    1: 'H',
    2: 'He',
    3: 'Li',
    4: 'Be',
    5: 'B',
    6: 'C',
    7: 'N',
    8: 'O',
    9: 'F',
    10: 'Ne',
    17: 'Cl'
}

class InvalidFileFormat(Exception):
    pass

class Atom:
    def __init__(self, atomic_number, coordinates):
        self.atomic_number = int(atomic_number)
        self.symbol = ELEMENTS[self.atomic_number]
        self.coordinates = np.array(coordinates, dtype=float)

    def to_jsonable(self):
        return {
            'atomic_number': self.atomic_number,
            'coordinates': self.coordinates.tolist()
        }

    def to_line(self):
        res = ' ' + ELEMENTS.get(self.atomic_number)
        np.set_printoptions(suppress=True)
        res += ' ' + np.array2string(self.coordinates)[1:-1]
        return res + '\n'

    def get_distance(self, other):
        return np.linalg.norm(other.coordinates - self.coordinates)
    
    # Returns an angle OTHER1—SELF—OTHER2
    def get_angle(self, other1, other2):
        v1 = other1.coordinates - self.coordinates
        v2 = other2.coordinates - self.coordinates
        cosine = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        return np.degrees(np.arccos(cosine))

    # Return a torsion angle between a1—a2 and a4—a3
    # Algorithm from https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python
    @staticmethod
    def get_torsion_angle(a1, a2, a3, a4):
        v1 = a1.coordinates - a2.coordinates
        v2 = a3.coordinates - a2.coordinates
        v3 = a4.coordinates - a3.coordinates
        
        v2 /= np.linalg.norm(v2)
        v = v1 - np.dot(v1, v2) * v2
        w = v3 - np.dot(v3, v2) * v2

        x = np.dot(v, w)
        y = np.dot(np.cross(v2, v), w)
        return np.degrees(np.arctan2(y, x))

class Molecule:
    def __init__(self, jsonDict=None):
        if jsonDict is None:
            self.atoms = []
        else:
            self.atoms = [Atom(i['atomic_number'], i['coordinates']) for i in jsonDict]

    def to_jsonable(self):
         return [atom.to_jsonable() for atom in self.atoms]

    def add_atom(self, atom):
        self.atoms.append(atom)

    def to_z_matrix(self, special={}):
        if len(self.atoms) == 0:
            return ''

        first = self.atoms[0]
        z_matrix = ' %s\n' % first.symbol
        if len(self.atoms) == 1:
            return z_matrix

        second = self.atoms[1]
        z_matrix += ' %s, 1, %f\n' % (second.symbol, second.get_distance(first))
        if len(self.atoms) == 2:
            return z_matrix

        third = self.atoms[2]
        z_matrix += ' %s, 1, %f, 2, %f\n' % (
            third.symbol,
            third.get_distance(first),
            first.get_angle(third, second)
        )
        for i in range(3, len(self.atoms)):
            a = self.atoms[i]
            if (i + 1) in special:
                string = ' %s, {}, %f, {}, %f, {}, %f\n'.format(*special[i + 1])
                ref = [self.atoms[k - 1] for k in special[i + 1]]
            else:
                string = ' %s, 1, %f, 2, %f, 3, %f\n'
                ref = [first, second, third]
            z_matrix += string % (
                a.symbol,
                a.get_distance(ref[0]),
                ref[0].get_angle(a, ref[1]),
                a.get_torsion_angle(a, ref[0], ref[1], ref[2])
            )
        return z_matrix

def find(f, search_for, shift=0, start_from=1, compare=lambda search_for, read_line: search_for == read_line):
    f.seek(shift, start_from)
    line = f.readline()
    while line and (not compare(search_for, line)):
        line = f.readline()
    if not line:
        return -1
    return f.seek(0, 1)

def split_line(line):
    return [el.replace('\n', '') for el in line.split(' ') if el.replace('\n', '')]

def skip_optimization_info(outfile):
    return skip_till_line('-- Stationary point found', outfile)

def skip_till_line(till, outfile):
    line = outfile.readline()
    while line:
        if till in line:
            return True
        line = outfile.readline()
    return False

def skip_n_lines(n, outfile):
    while n > 0:
        n -= 1
        outfile.readline()

def read_optimized_coords(outfile):
    raw_coords = []
    line = outfile.readline()
    while line and (LINE_SEPARATOR not in line):
        raw_coords.append(line)
        line = outfile.readline()
    return raw_coords

def get_optimized_coords(outfile, ignore_error=False):
    ST_OR = 'Standard orientation:'
    f = open(outfile)
    if not ignore_error:
        if not skip_optimization_info(f):
            f.close()
            raise InvalidFileFormat('No optimization info')
    pos = find(f, ST_OR, compare=lambda search_for, read_line: search_for in read_line)
    if pos == -1:
        f.close()
        raise InvalidFileFormat('No coordinates data')
    while pos != -1:
        prev_pos = pos
        pos = find(f, ST_OR, compare=lambda search_for, read_line: search_for in read_line)
    f.seek(prev_pos, 0)
    skip_n_lines(4, f)
    raw_coords = read_optimized_coords(f)
    molecule = Molecule()
    for line in raw_coords:
        parts = split_line(line)
        atom = Atom(parts[1], parts[3:6])
        molecule.add_atom(atom)
    f.close()
    return molecule

def get_free_gibbs_energy(outfile):
    reached = False
    f = open(outfile)
    for line in f:
        if 'Sum of electronic and thermal Free Energies' in line:
            reached = True
            break
    f.close()
    if not reached:
        raise InvalidFileFormat('No information on thermal properties')
    return float(split_line(line)[-1])

def get_scan_results(outfile):
    f = open(outfile)
    reached = skip_till_line('Summary of the potential surface scan:', f)
    if not reached:
        raise InvalidFileFormat('No scan results info')
    headers = split_line(f.readline())
    skip_n_lines(1, f)
    data = [[] for i in headers]
    for line in f:
        if '----' in line:
            break
        d = split_line(line)
        for i in range(len(headers)):
            data[i].append(d[i])
    result = {
        headers[0]: np.array(data[0], dtype=int)
    }
    for i in range(1, len(headers)):
        result[headers[i]] = np.array(data[i], dtype=float)
    return result

def write_coords_file(molecule, file, charge=0, multiplicity=1):
    f = open(file, mode='w', encoding='ascii')
    f.write(str(charge) + ' ' + str(multiplicity) + '\n')
    for atom in molecule.atoms:
        f.write(atom.to_line())
    f.close()