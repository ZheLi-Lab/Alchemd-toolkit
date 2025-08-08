import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
from abc import ABC, abstractmethod 
from parse import parse 
import sys 
from optparse import OptionParser 
from scipy import integrate, interpolate, signal 
from scipy.interpolate import CubicSpline 
import re 
import os 
from sklearn.linear_model import Ridge, Lasso, RidgeCV, LassoCV 
from sklearn.pipeline import make_pipeline 
from sklearn.preprocessing import StandardScaler 
import matplotlib 
matplotlib.use('agg') 

class optParser():  
    def __init__(self, fakeArgs): 
        parser = OptionParser() 
        parser.add_option('-f', '--free_ene_csv', dest='free_ene_csv', help='The name of the free energy csv file to be analyzed', default='FEPfree_ene.csv') 
        parser.add_option('-r','--lig_refer_csv', dest='lig_refer_csv', help='The reference file (generally ligand free_ene file) to generate ridge model.',default=None) 
        parser.add_option('-l', '--lambda_col_name', dest='lambda_col_name', help='The name of the column that record lambda.', default='delta_A_what_to_what') 
        parser.add_option('-e', '--fe_col_name', dest='fe_col_name', help='The name of the column that record free energy.', default='free_energy(kcal/mol)') 
        parser.add_option('-p', '--ifsavepng', dest='ifsavepng', help='Whether to output png file.', default=True) 
        parser.add_option('-m', '--fe_mode', dest='fe_mode', help='The mode of the free energy csv file. BAR or FEP', default='FEP') 
        parser.add_option('-w', '--cal_wins', dest='cal_wins', help='The filename that give the windows information for calculation. \ If you give the string of "all", the program will not filter any data.', default='all') 
        parser.add_option('-z', '--md_pack', dest='md_pack', help='Md simulation pack, openmm or amber.', default='openmm') 
        if fakeArgs: 
            self.option, self.args = parser.parse_args(fakeArgs) 
        else: 
            self.option, self.args = parser.parse_args() 

class dG_dlambda(): 
    def __init__(self, dG_data, md_pack='openmm', fe_mode='FEP'): 
        self.fe_mode = fe_mode 
        self.df_csv = dG_data
        self.md_pack = md_pack

    @staticmethod
    def cal_delta_lambda(lambda_value_1, lambda_value_2):
        lambda_value_1 = np.array(lambda_value_1)
        lambda_value_2 = np.array(lambda_value_2)
        dd_coords = abs(np.subtract(lambda_value_1,lambda_value_2))
        ABFE_labels = ['restraint', 'charge', 'vdw']
        if type(dd_coords)==np.float64:
            label = 'amber_RBFE'
            delta_lambda = np.around(dd_coords, decimals=3) 
            lambda_info_x = np.around((lambda_value_1+lambda_value_2)/2, decimals=4)
        else:
            nonzero_indices = np.argwhere(dd_coords != 0).flatten()
            nonzero_count = len(nonzero_indices)
            if nonzero_count==1:
                delta_lambda = np.around(dd_coords[nonzero_indices[0]], decimals=3)
                if nonzero_indices[0] < len(ABFE_labels):
                    label = ABFE_labels[nonzero_indices[0]]
                    lambda_info_x = np.around((lambda_value_1[nonzero_indices[0]]+lambda_value_2[nonzero_indices[0]])/2, decimals=4)
                else:
                    label = 'unknown'
                    lambda_info_x = None
            elif nonzero_count==2 and lambda_value_1[1]==lambda_value_1[2] and lambda_value_2[1]==lambda_value_2[2]:
                label = 'openmm_RBFE'
                delta_lambda = np.around(dd_coords[nonzero_indices[0]], decimals=3)
                lambda_info_x = np.around((lambda_value_1[nonzero_indices[0]]+lambda_value_2[nonzero_indices[0]])/2, decimals=4)
            else:
                label = 'mix'
                delta_lambda = np.around(dd_coords.sum(), decimals=3) 
                lambda_info_x = None
        return delta_lambda, label, lambda_info_x
    
    @staticmethod
    def get_lambda_start_end_coords(delta_A_what_to_what, md_pack):
        pattern = '({},{},{}) to ({},{},{})'
        result = parse(pattern, delta_A_what_to_what)
        if result:
            start_coords = (float(result[0]), float(result[1]), float(result[2]))
            end_coords = (float(result[3]), float(result[4]), float(result[5]))
        else:
            pattern = '{} to {}'
            result = parse(pattern, delta_A_what_to_what)
            if result:
                start_coords = float(result[0])
                end_coords = float(result[1])
            else:
                raise ValueError("delta_A_what_to_what can't match pattern '({},{},{}) to ({},{},{})' or '{} to {}'. Please check the format.")
        return start_coords, end_coords

    def get_dG_dlambda(self, lambda_col_name, free_ene_col_name):
        start, end = type(self).get_lambda_start_end_coords(self.df_csv[lambda_col_name].values[-1],self.md_pack)
        lambda_col_exceptLastLine = str(self.df_csv[lambda_col_name].values[0:-1])
        if str(start) in lambda_col_exceptLastLine and str(end) in lambda_col_exceptLastLine:
            self.df_csv = self.df_csv.iloc[0:len(self.df_csv.index)-1,:]
        # print(self.df_csv)
        dG_dlambda_values = []
        labels = []
        lambda_infos = []
        for index, row in self.df_csv.iterrows():
            lambda_col_value = row[lambda_col_name]
            free_ene_col_value = row[free_ene_col_name]
            start_coords, end_coords = type(self).get_lambda_start_end_coords(lambda_col_value, self.md_pack)
            dlambda, label, lambda_info_x = type(self).cal_delta_lambda(np.array(start_coords), np.array(end_coords))
            dG_dlambda = free_ene_col_value/dlambda
            dG_dlambda_values.append(dG_dlambda)
            labels.append(label)
            lambda_infos.append(lambda_info_x)
        self.df_csv = self.df_csv.assign(dG_dlambda=dG_dlambda_values)
        self.df_csv = self.df_csv.assign(lambda_labels=labels)
        self.df_csv = self.df_csv.assign(lambda_info=lambda_infos)

    def plot_dG_dlambda(self, df, postfix=None):
        '''
        plotting dG/dlambda individually by different 'lambda_labels' in given df.
        '''
        grouped_df = df.groupby('lambda_labels', sort=False)
        for key, value in grouped_df:
            if 'RBFE' in key or key=='charge' or key=='vdw':
                if postfix is not None:
                    pngname = key + postfix
                    self.plotting(value, pngname)

    def plotting(self, df, png_file=None):
        plt.clf()
        # Set bar chart colors
        color_dict = {
            'restraint': 'tab:red',
            'charge': 'tab:blue',
            'vdw': 'tab:green',
            'mix': 'tab:orange',
            'amber_RBFE': '#EDB120',
            'openmm_RBFE': '#7E2F8E'
        }
        # plt.figure(figsize=(10,10))
        # Draw bar chart
        plt.bar(df.index, df['dG_dlambda'], color=[color_dict[l] for l in df['lambda_labels']])
        # Draw line
        plt.plot(df.index, df['dG_dlambda'], color='black')
        # Set x-axis label
        plt.xticks(df.index, df['delta_A_what_to_what'], rotation=90)
        plt.ylabel(r'$\delta G_{cal}/\delta \lambda $ (kcal/mol)', )#y-axis label
        plt.xlabel(r'$\lambda$ info',)#x axis label
        
        # plt.minorticks_on()#Enable minor ticks
        plt.tick_params(which='major',width=3, length=6)#Set major tick size
        # plt.tick_params(which='minor',width=2, length=4)#Set minor tick size
        ax=plt.gca();#Get axis handle
        ax.spines['bottom'].set_linewidth(4);#Set bottom axis thickness
        ax.spines['left'].set_linewidth(4);#Set left axis thickness
        ax.spines['right'].set_linewidth(4);#Set right axis thickness
        ax.spines['top'].set_linewidth(4);#Set top axis thickness
        plt.tight_layout()
        if png_file is not None:
            plt.savefig(f'{png_file}.png', format="png",dpi=600, transparent=True)
        else:
            plt.show()
        plt.clf()

class DATA():
    def __init__(self, data, cal_lst):
        self.cal_lst = cal_lst
        if self.cal_lst != 'all':
            # print(f'cal_win_lst is:{self.cal_lst}')
            data = self.extract_part(data, self.cal_lst)
            # print(f'data:{data}')
        self.df = data.drop_duplicates(['lambda_info']).dropna(subset=['dG_dlambda'])
        self.lambda_info = np.array(self.df['lambda_info']) 
        self.dG_dlambda = np.array(self.df['dG_dlambda'])

        self.lambda_label = list(self.df['lambda_labels'])[0]
        self.png_prefix = 'CSpline-'+self.lambda_label
    def extract_part(self, df, cal_lst):
        for i in range(len(cal_lst)):
            cal_lst[i] = str(cal_lst[i])
        df = df.loc[df['lambda_value'].isin(cal_lst)]
        return df

class GetCurve():
    def __init__(self, data, cal_lst, method, n=3):
        self.data_obj = DATA(data, cal_lst)
        self.df = self.data_obj.df
        ###x: 1-D array of lambda values, and should be in strictly increasing order.
        self.x = self.data_obj.lambda_info
        ###y: An array with arbitrary number of dimensions, but every axis should have the same length with x.
        self.y = self.data_obj.dG_dlambda
        self.lambda_label = list(self.df['lambda_labels'])[0]
        self.method = method
        if self.method == 'cspline':
            self.curve_func = self.get_cubic_spline(self.x, self.y)
        elif self.method == 'polynomial':
            self.curve_func = self.get_polyfunc(self.x, self.y, n)
        else:
            print(f'Error: {method} is not supported, please choose from [\'cspline\', \'polynomial\'].')
        # self.cspline = self.get_cubic_spline(self.x, self.y)
        # self.polyfunc = self.get_polyfunc(self.x, self.y, n)


    def get_cubic_spline(self, x, y):
        '''
        cspline: object, like a function, y_predict = cspline(x_toPredict).
        '''
        cspline = CubicSpline(x, y, bc_type=((2, 1), (2, 1)))
        return cspline

    def cal_cubic_inte(self, cspline, ifplot=False, postfix='dGdlambda', dot_num=100):
        '''
        Return
        ----------
        inte: The integrate of CubicSpline generated by value of lambda_info and dG_dlambda.
        '''
        cubic_inte = cspline.integrate(0,1)
        # print(f'integrate of cubic spline is:{inte}')
        if ifplot:
            self.plot_curve(cspline, postfix, dot_num)
        return float(cubic_inte)
    
    def get_polyfunc(self, x, y, n=3):
        #polyfit
        coefficients = np.polyfit(x, y, n, full=False, w=None, cov=False)
        polyfunc = np.poly1d(coefficients)
        # print(self.polyfunc)
        return polyfunc

    def cal_polyfunc_inte(self, polyfunc, ifplot=False, postfix='dGdlambda', dot_num=100):
        intefunc = np.polyint(polyfunc, m=1)
        # print(f'integrate is: {intefunc(1)}')
        poly_inte = intefunc(1)-intefunc(0)
        # print(f'integrate of polynomial is:{inte}')
        if ifplot:
            self.plot_curve(polyfunc, postfix, dot_num)
        return float(poly_inte)

    def cal_y_predict(self, curve_func, x):
        '''
        calculate y_predict by given fitting function and input x.
        '''
        return curve_func(x)

    def plot_curve(self, curve_func, postfix='dGdlambda', dot_num=100):
        x2 = np.linspace(0,1,dot_num).round(2)
        y2 = self.cal_y_predict(curve_func, x2)
        ##plot original data and fitting curve
        plt.clf()
        plt.plot(self.x, self.y, 'o', color='black', markersize=4, label='data')
        plt.plot(x2, y2, color='black', markersize=3, label=self.method)
        plt.legend()
        ## Set axis labels
        plt.xticks(self.x, self.df['delta_A_what_to_what'], rotation=90)
        plt.ylabel(r'$\delta G_{cal}/\delta \lambda $ (kcal/mol)', )#y-axis label
        plt.xlabel(r'$\lambda$ info',)#x-axis label
        # plt.minorticks_on()#Enable minor ticks
        plt.tick_params(which='major',width=3, length=6)#Set major tick size
        # plt.tick_params(which='minor',width=2, length=4)#Set minor tick size
        ax=plt.gca();#Get axis handle
        ax.spines['bottom'].set_linewidth(4);#Set bottom axis thickness
        ax.spines['left'].set_linewidth(4);#Set left axis thickness
        ax.spines['right'].set_linewidth(4);#Set right axis thickness
        ax.spines['top'].set_linewidth(4);#Set top axis thickness
        plt.tight_layout()
        plt.savefig(f'{self.method}-{self.lambda_label}_{postfix}.png',dpi=600)
        plt.clf()
        return y2

class PolyFit():
    def __init__(self, model_data, to_fit_data):
        self.model_data = model_data ##the csv to generate polynomial model
        self.to_fit_data = to_fit_data ##the target csv to fit with modeled polynomial curve

    def gen_polyfit_model(self, x, y, n=3, ifplot=False, prefix='model'):
        x = np.array(x) ##lambda value as x axis
        y = np.array(y) ##dG/dlambda value as y axis
        #polyfit
        coefficients = np.polyfit(x, y, n, full=False, w=None, cov=False)
        polyfunc = np.poly1d(coefficients)
        print(polyfunc)
        if ifplot:
            #plot original data and polynomial curve
            plt.clf()
            y_new = np.polyval(polyfunc,x)
            plt.plot(x, y, '*', markersize=4)
            plt.plot(x, y_new, markersize=4)
            plt.savefig(f'{prefix}_polynomialCurve.png',dpi=600)
            plt.close()

class RidgeFit():
    def __init__(self, model_data, model_win_lst, to_fit_data, to_cal_lst, model_mode='cspline', n=3):
        '''
        model_mode: a string, determine the method to build model curve. {'cspline':cubic spline | 'polynomial':polynomials}
        '''
        self.model_mode = model_mode
        self.n = n
        self.to_fit_data_obj = DATA(to_fit_data, to_cal_lst)
        self.to_fit_x = self.to_fit_data_obj.lambda_info
        self.to_fit_y = self.to_fit_data_obj.dG_dlambda
        self.to_fit_X = self.gen_X(self.to_fit_x, model_data, model_win_lst)
        self.ridge_model = self.get_ridge_model(self.to_fit_X, self.to_fit_y)
        
    def gen_X(self, x, model_data, model_win_lst):
        '''
        Generate the array X used for Ridge regression
        Return
        ----------
        X: array, [[1, x1, f(x1)_predictByModel],[1, x2, f(x2)_predictByModel], ...,[1, xn, f(xn)_predictByModel]].
        '''
        model_curve_obj = GetCurve(model_data, model_win_lst, self.model_mode, self.n)
        self.model_data_obj = model_curve_obj.data_obj
        predict_y = model_curve_obj.curve_func(x)
        X = np.array([[1]*len(x), x, predict_y]).transpose()
        return X

    def get_ridge_model(self, X, y):
        ###calculate the best alpha by cross validation.
        model_alphas = np.logspace(-3, 2, 20)
        CVestimator = RidgeCV(alphas=model_alphas)
        # ridge_best_Lambda = self.cal_best_Lambda(X, y, CVestimator)
        ridge_best_Lambda = 0.01
        ridge_estimator = Ridge(alpha=ridge_best_Lambda)
        ridge_model = make_pipeline(StandardScaler(), ridge_estimator)
        ridge_model.fit(X, y)
        return ridge_model
        ## no normalization
        # ridge_estimator.fit(self.X, self.y)
        # return ridge_estimator

    def cal_best_Lambda(self, X, y, CVestimator):   
        esti_type = re.match(".{7}",f'{CVestimator}').group(0)
        CVestimator.fit(X, y)
        best_Lambda = CVestimator.alpha_
        print(f'{esti_type}_best_Lambda is:{best_Lambda}')
        return best_Lambda
    
    def cal_y_predict(self, ridge_model, x):
        X = self.gen_X(x, self.model_data_obj.df, self.model_data_obj.cal_lst)
        y_pred = ridge_model.predict(X) ##array [[y1],[y2],[y3],...,[yn]]
        return y_pred.T ##array [y1, y2, y3, ..., yn]

    def cal_ridge_inte(self, ridge_model, ifplot=False, postfix='dGdlambda', dot_num=100):
        pred_x = np.linspace(0,1,dot_num).round(2)
        y_new = self.cal_y_predict(ridge_model, pred_x)
        ridge_inte = integrate.trapz(y_new, pred_x)
        if ifplot:
            self.plot_ridge(pred_x, y_new, postfix)
        return ridge_inte
    
    def plot_ridge(self, pred_x, y_predict, postfix='dGdlambda'):
        plt.clf()
        plt.plot(self.to_fit_x, self.to_fit_y, 'o', color='black', markersize=4, label='data')
        plt.plot(pred_x, y_predict, color='black', markersize=3, label='ridge')
        plt.legend()
        ## Set axis labels
        plt.xticks(self.to_fit_x, self.to_fit_data_obj.df['delta_A_what_to_what'], rotation=90)
        plt.ylabel(r'$\delta G_{cal}/\delta \lambda $ (kcal/mol)', )#y-axis label
        plt.xlabel(r'$\lambda$ info',)#x-axis label
        # plt.minorticks_on()#Enable minor ticks
        plt.tick_params(which='major',width=3, length=6)#Set major tick size
        # plt.tick_params(which='minor',width=2, length=4)#Set minor tick size
        ax=plt.gca();#Get axis handle
        ax.spines['bottom'].set_linewidth(4);#Set bottom axis thickness
        ax.spines['left'].set_linewidth(4);#Set left axis thickness
        ax.spines['right'].set_linewidth(4);#Set right axis thickness
        ax.spines['top'].set_linewidth(4);#Set top axis thickness
        plt.tight_layout()
        plt.savefig(f'ridge_{self.model_mode}-{self.to_fit_data_obj.lambda_label}_{postfix}.png',dpi=600)
        plt.clf()
        return y_predict

class ANA():
    def __init__(self, get_csv_obj, fe_mode, ifplot=False):
        self.get_csv_obj = get_csv_obj
        self.fe_mode = fe_mode
        self.ifplot = ifplot
        if self.ifplot:
            self.get_csv_obj.plot_dG_dlambda(get_csv_obj.df_csv,f'{fe_mode}_dGdlambda')

    
    def ana_fitting(self, to_fit_df, cal_wins, model_data=None, method='cspline', n=3, index_name='cal_win_all', dot_num=100):
        ## model data is the ligand_vdw_dGdlambda for ABFE or ligand_dGdlambda for RBFE.
        ## to_fit_data is complex_vdw_dGdlambda DataFrame for ABFE or complex_dGdlambda DataFrame for RBFE.
        if model_data:
            ridge_obj = RidgeFit(model_data, 'all', to_fit_df, cal_wins, method, n)
            inte = ridge_obj.cal_ridge_inte(ridge_obj.ridge_model, self.ifplot, f'{self.fe_mode}{index_name}',dot_num)
        else:
            curve_obj = GetCurve(to_fit_df, cal_wins, method, n)
            if method == 'cspline':
                inte = curve_obj.cal_cubic_inte(curve_obj.curve_func, self.ifplot, f'{self.fe_mode}_{index_name}',dot_num)
            elif method == 'polynomial':
                inte = curve_obj.cal_polyfunc_inte(curve_obj.curve_func, self.ifplot, f'{self.fe_mode}_{index_name}',dot_num)
            else:
                print(f'Error: {method} is not supported, please choose from [\'cspline\', \'polynomial\'].')
        return inte

    def cal_inte(self, get_csv_obj, cal_wins, method, n, index_name, ligand_file=None, dot_num=100):
        df_csv = get_csv_obj.df_csv
        print(df_csv)
        grouped_df = df_csv.groupby('lambda_labels', sort=False)
        groups_index = grouped_df.size().index
        inte_result = []
        if len(groups_index)==1:
            try:
                ligand_data = pd.read_csv(ligand_file, sep='|')
            except:
                ligand_data = None
            inte = self.ana_fitting(df_csv, cal_wins, ligand_data, method, n, index_name, dot_num)
            inte_result.append(inte)
        elif 'vdw' in groups_index or 'charge' in groups_index:
            if 'charge' in groups_index:
                df_crg = grouped_df.get_group('charge')
                inte_crg = self.ana_fitting(df_crg, cal_wins, ligand_file, method, n, index_name, dot_num)
                inte_result.append(inte_crg)
                df_vdw = grouped_df.get_group('vdw')
            if 'vdw' in groups_index:
                try:
                    ligand_data = pd.read_csv(ligand_file, sep='|').groupby('lambda_labels', sort=False).get_group('vdw')
                except:
                    ligand_data = None
                inte_vdw = self.ana_fitting(df_vdw, cal_wins, ligand_data, method, n, index_name, dot_num)
                inte_result.append(inte_vdw)
        else:
            print('Lambda labels are morethan one, but charge and vdw data is not found. Please check your data if need.')
        return np.array(inte_result)

    def do_ana(self, cal_wins, fe_col_name, method='cspline', n=3, ligand_file=None, dot_num=100):
        if 'restraint' in self.get_csv_obj.df_csv['lambda_labels'].values:
            dG_restraint = float(self.get_csv_obj.df_csv[fe_col_name][np.where(self.get_csv_obj.df_csv['lambda_labels']=='restraint')[0]])
        else:
            dG_restraint = 0

        df_result = pd.DataFrame()
        inte_result = self.cal_inte(self.get_csv_obj, 'all', method, n, 'cal_win_all', ligand_file, dot_num)
        single_df = pd.DataFrame(data=[[dG_restraint, inte_result, dG_restraint+inte_result.sum()]], \
                                 columns=['restr', 'integrate_result', 'sum'], index=['cal_win_all'])
        df_result = pd.concat([df_result, single_df], axis=0)        
        if cal_wins != 'all':
            for idx_ in range(0, len(cal_wins)):
                win_lst = cal_wins[idx_]
                inte_result = self.cal_inte(self.get_csv_obj, win_lst, method, n, f'cal_win_{idx_}', ligand_file, dot_num)
                single_df = pd.DataFrame(data=[[dG_restraint, inte_result, dG_restraint+inte_result.sum()]],\
                                         columns=['restr', 'integrate_result', 'sum'], index=[f'cal_win_{idx_}'])
                df_result = pd.concat([df_result, single_df], axis=0)
        df_result.index.names=['cal_win']
        df_result.to_csv(f'fitting_{self.fe_mode}.csv')

if __name__ == '__main__':
    opts = optParser('') 
    free_ene_csv = str(opts.option.free_ene_csv)
    lambda_col_name = str(opts.option.lambda_col_name)
    fe_col_name = str(opts.option.fe_col_name)
    fe_mode = str(opts.option.fe_mode)
    ifsavepng = opts.option.ifsavepng
    ligand_refer_file = opts.option.lig_refer_csv
    md_pack = str(opts.option.md_pack)

    ## get cal_win_lst from the specific file.
    if opts.option.cal_wins == 'all':
        cal_win_lst='all'
    else:
        with open(opts.option.cal_wins, 'r') as f:
            f_content = f.readlines()
            cal_win_lst = [eval(line.strip()) for line in f_content]
            if not isinstance(cal_win_lst, list):
                print('Error: The content of the win_lst file is not a list!')
                sys.exit()
            print(len(cal_win_lst))

    dG_data = pd.read_csv(free_ene_csv, sep='|')
    get_csv_obj = dG_dlambda(dG_data, md_pack, fe_mode)
    get_csv_obj.get_dG_dlambda(lambda_col_name, fe_col_name)
    ts = ANA(get_csv_obj, fe_mode, ifsavepng)
    ts.do_ana(cal_win_lst, fe_col_name,  'cspline', 3, ligand_refer_file, 100)
