import os
import re
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
plt.rcParams['figure.figsize'] = 12,5

# Read in command line arguement
if len(sys.argv) > 1:
    to_plot = sys.argv[1]
else:
    to_plot = False

print(os.listdir('results'))

files = os.listdir('results')
results = [r for r in files if r.startswith('r')]

catch = r'results\_([0-9]*)\_([0-9]*)\_([beam0-9]*)\_(0\.[0-9]*)\.csv'

groups = [ re.match( catch,r ).groups() for r in results ]
for g in groups:
    print(g)

n_vals = [ int( g[0] )  for g in groups ]
p = groups[0][1]

if str(groups[0][2]) == "beam":
    beam = True
    beam_h = [g[-1] for g in groups]
else:
    beam = False

# Define n, h, and f values
h_vals = [1./n for n in n_vals]
f_ords = [0,1,2]
f_labs = [r'$f(x) = c$',r'$f(x)=x$',r'$f(x)=x^2$']

# Load results from program
errors = pd.read_csv('results/errors.csv')
error_str = r'$||u-u^h||_{L^2}$'

# Plot the approximate solutions
if to_plot:
    if beam:
        groups = errors.groupby(['h','n']) 
        for i,h in enumerate(beam_h)
        for i,n in enumerate(n_vals):
            # Get error from dataframe
            error = groups.get_group((f_ord, n))['error'].values[0]
            error = round(error,10)
            fname = 'results/results_' + str(n) + '_' + str(f_ord) + '.csv'

            # Load solutions for this choice of n and f
            res = pd.read_csv(fname, index_col = 0)
            d = res.index
            uh = res['uh']
            u = res['u']

            # Plot the solutions
            plt.subplot(1,4,i+1)        
            if i==0: plt.ylabel(r'$u^h(x)$')
            plt.plot(d, u, 'C0', label=r'$u(x)$',linewidth=2.5)
            plt.plot(d, uh, ':k', label=r'$u^h(x)$',linewidth=2.5)

            # Fomrate the plot
            plt.title('n = {}\n{} = {}'.format(n, error_str, error))
            plt.legend()
            plt.xlabel('x')
            plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    for f_ord, f_lab in zip(f_ords,f_labs):
        for i,n in enumerate(n_vals):
            # Get error from dataframe
            groups = errors.groupby(['forder','n']) 
            error = groups.get_group((f_ord, n))['error'].values[0]
            error = round(error,10)
            fname = 'results/results_' + str(n) + '_' + str(f_ord) + '.csv'

            # Load solutions for this choice of n and f
            res = pd.read_csv(fname, index_col = 0)
            d = res.index
            uh = res['uh']
            u = res['u']

            # Plot the solutions
            plt.subplot(1,4,i+1)        
            if i==0: plt.ylabel(r'$u^h(x)$')
            plt.plot(d, u, 'C0', label=r'$u(x)$',linewidth=2.5)
            plt.plot(d, uh, ':k', label=r'$u^h(x)$',linewidth=2.5)

            # Fomrate the plot
            plt.title('n = {}\n{} = {}'.format(n, error_str, error))
            plt.legend()
            plt.xlabel('x')
            plt.tight_layout(rect=[0, 0.03, 1, 0.95])

        # Display the plot
        plt.suptitle(f_lab)
        plt.savefig('figures/solutions_f_order_{}'.format(f_ord)+'.png')
        plt.show()


# Plot log-log error plots
for i in range(len(f_ords)):
    plt.subplot(1,3,i+1)

    # Get errors for this order of f
    error_vals = errors[errors['forder']==f_ords[i]]['error'].values        

    # Plot the errors against h
    plt.loglog(h_vals,error_vals)

    plt.loglog(h_vals,np.array(h_vals)**2,label=r'$h^2$')
    plt.legend()

    # calculate slope of error plot. 
    rise = np.log(error_vals[-1]) - np.log(error_vals[0])
    run = np.log(h_vals[-1]) - np.log(h_vals[0])
    slope = round( rise / run , 5 )

    # Format the plot
    plt.title('slope: {}\n'.format(slope) + f_labs[i])
    plt.xlabel('h')
    if i==0: plt.ylabel(error_str)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    

# Display the plot
plt.suptitle('Log-Log Error Plot')
plt.savefig('figures/log_log_errors.png')
plt.show()

