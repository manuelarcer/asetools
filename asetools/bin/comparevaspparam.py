#!/usr/bin/env python

def read_outcar(filename):
    """Reads an OUTCAR file and extracts parameters into a dictionary."""
    parameters = {}
    keywords = ["POTCAR", "A1", "A2", "A3", "NKPTS", "NKDIM", "NBANDS", "NEDOS", "NIONS", "LDIM", "LMDIM", "NPLWV",
                "IRMAX", "IRDMAX", "NGX", "NGY", "NGZ", "NGXF", "NGYF", "NGZF", "NWRITE", "PREC", "ISTART", "ICHARG",
                "ISPIN", "LNONCOLLINEAR", "LSORBIT", "INIWAV", "LASPH", "METAGGA", "ENCUT", "ENINI", "ENAUG", "NELM",
                "NELMIN", "NELMDL", "EDIFF", "LREAL", "NLSPLINE", "LCOMPAT", "GGA_COMPAT", "LMAXPAW", "LMAXMIX",
                "VOSKOWN", "ROPT", "EDIFFG", "NSW", "NBLOCK", "KBLOCK", "IBRION", "NFREE", "ISIF", "IWAVPR", "ISYM",
                "LCORR", "POTIM", "TEIN", "TEBEG", "TEEND", "SMASS", "NPACO", "APACO", "PSTRESS", "IALGO", "LDIAG",
                "LSUBROT", "TURBO", "IRESTART", "NREBOOT", "NMIN", "EREF", "IMIX", "AMIX", "BMIX", "AMIX_MAG",
                "BMIX_MAG", "AMIN", "WC", "INIMIX", "MIXPRE", "MAXMIX", "LMONO", "LDIPOL", "IDIPOL", "EPSILON",
                "GGA", "LEXCH", "LHFCALC", "LHFONE", "AEXX", "LEPSILON", "LRPA", "LNABLA", "LVEL", "LINTERFAST",
                "KINTER", "CSHIFT", "OMEGAMAX", "DEG_THRESHOLD", "RTIME", "DFIELD", "ORBITALMAG", "LCHIMAG", "DQ",
                "LLRAUG"]
    with open(filename, 'r') as file:
        for line in file:
            # Extract parameters based on known keywords
            for keyword in keywords:
                if line.startswith(keyword):
                    parameters[keyword] = line.strip()
    return parameters

def compare_parameters(params1, params2):
    """Compares two dictionaries of parameters and prints differences."""
    keys = set(params1.keys()).union(params2.keys())
    for key in keys:
        if params1.get(key) != params2.get(key):
            print(f"Difference in {key}:")
            print(f"File 1: {params1.get(key)}")
            print(f"File 2: {params2.get(key)}\n")

def main():
    file1 = input("Enter the path to the first OUTCAR file: ")
    file2 = input("Enter the path to the second OUTCAR file: ")
    params1 = read_outcar(file1)
    params2 = read_outcar(file2)
    compare_parameters(params1, params2)

if __name__ == "__main__":
    main()
