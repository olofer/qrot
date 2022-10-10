import numpy as np 
import argparse 
import matplotlib.pyplot as plt

def load_csv_and_colnames(csvname):
  with open(csvname) as f:
    first_line = f.readline().strip()
  Xnames = [e.strip() for e in first_line.split(',')]
  X = np.loadtxt(csvname, skiprows = 1, delimiter = ',')
  assert X.shape[1] == len(Xnames)
  return X, Xnames

if __name__ == '__main__':
  parser = argparse.ArgumentParser()

  parser.add_argument('--trace-files', type = str, nargs = '+', help = 'name of log files (CSV) to load')
  parser.add_argument('--figext', type = str, default = 'pdf', help = 'figure file extension (e.g. pdf or png)')
  parser.add_argument('--dpi', type = float, default = 250.0, help = 'print to file PNG option')

  args = parser.parse_args()

  assert len(args.trace_files) > 0
  print(args.trace_files)

  table = {}
  for f in args.trace_files:
    print(f)
    X, Xnames = load_csv_and_colnames(f)
    table[f] = (Xnames, X)

  # spin stability overview plot
  for f in args.trace_files:
    Xnames, X = table[f]
    plt.plot(X[:, Xnames.index('time')], X[:, Xnames.index('costheta')], 
             linewidth = 2, alpha = 0.50, label = f.split('.')[0])
  plt.xlabel('time [sec]')
  plt.ylabel('$\\cos(\\theta)$')
  plt.grid(True)
  plt.title('Stability of spins')
  plt.legend()
  plt.tight_layout()
  plt.savefig('qrot-plot-theta.{}'.format(args.figext), dpi = args.dpi)
  plt.close()

  # conservation of energy
  for f in args.trace_files:
    Xnames, X = table[f]
    plt.plot(X[:, Xnames.index('time')], X[:, Xnames.index('K')], 
             linewidth = 2, alpha = 0.50, label = f.split('.')[0])
  plt.xlabel('time [sec]')
  plt.ylabel('Kinetic energy [J]')
  plt.grid(True)
  plt.title('Conservation of energy')
  plt.legend()
  plt.tight_layout()
  plt.savefig('qrot-plot-K.{}'.format(args.figext), dpi = args.dpi)
  plt.close()

  # angular velocity components
  for f in args.trace_files:
    Xnames, X = table[f]
    for o in ['omegax', 'omegay', 'omegaz']:
      plt.plot(X[:, Xnames.index('time')], X[:, Xnames.index(o)], 
               linewidth = 2, alpha = 0.50, label = o)
    plt.xlabel('time [sec]')
    plt.ylabel('[rad/s]')
    plt.grid(True)
    plt.title('Angular velocity components')
    plt.legend()
    plt.tight_layout()
    name_wo_ext = f.split('.')[0]
    plt.savefig('qrot-plot-omega-{}.{}'.format(name_wo_ext, args.figext), dpi = args.dpi)
    plt.close()

  # angular momentum components
  for f in args.trace_files:
    Xnames, X = table[f]
    for o in ['Lx', 'Ly', 'Lz']:
      plt.plot(X[:, Xnames.index('time')], X[:, Xnames.index(o)], 
               linewidth = 2, alpha = 0.50, label = o)
    plt.xlabel('time [sec]')
    plt.ylabel('[kg*m^2/s]')
    plt.grid(True)
    plt.title('Angular momentum components')
    plt.legend()
    plt.tight_layout()
    name_wo_ext = f.split('.')[0]
    plt.savefig('qrot-plot-L-{}.{}'.format(name_wo_ext, args.figext), dpi = args.dpi)
    plt.close()
