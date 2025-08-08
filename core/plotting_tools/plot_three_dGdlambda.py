from dG_dlambda import *

lambda_col_name = "delta_A_what_to_what"
fe_col_name="free_energy(kcal/mol)"

class All_sys_data():
    def __init__(self, csv_dict):
        self.df_dict = {}
        for key, value in csv_dict.items():
            dG_dlamda_obj = ABFE_dG_dlambda(value, 'BAR')
            dG_dlamda_obj.get_dG_dlambda(lambda_col_name, fe_col_name)
            df_csv = dG_dlamda_obj.df_csv.drop(dG_dlamda_obj.df_csv[dG_dlamda_obj.df_csv['lambda_labels']=='mix'].index)
            df_crg = df_csv[df_csv['lambda_labels']=='charge']
            self.df_dict[key] = df_crg
    
    def plot_dG_dlambda(self,png_file):
        plt.clf()
        for key in self.df_dict.keys():
            df = self.df_dict[key]
            plt.plot(df.index, df['dG_dlambda'], label=key)
        # plt.xticks(df.index, [ np.around(i, 2) for i in np.arange(0, 1.0, 0.01)], rotation=90, fontsize=5)
        plt.ylabel(r'$\delta G_{cal}/\delta \lambda $ (kcal/mol)', )#y-axis label
        plt.xlabel(r'$\lambda$ info',)#x-axis label
        plt.legend() 
        # plt.minorticks_on()#Enable minor ticks
        # plt.tick_params(which='major',width=3, length=6)#Set major tick size
        # plt.tick_params(which='minor',width=2, length=4)#Set minor tick size
        ax=plt.gca();#Get axis handle
        ax.spines['bottom'].set_linewidth(4);#Set bottom axis thickness
        ax.spines['left'].set_linewidth(4);#Set left axis thickness
        ax.spines['right'].set_linewidth(4);#Set right axis thickness
        ax.spines['top'].set_linewidth(4);#Set top axis thickness
        plt.tight_layout()
            # Display figure
        if png_file is not None:
            plt.savefig(png_file, format="png",dpi=600, transparent=True)
        else:
            plt.show()

if __name__ == '__main__':
    csv_dict = {'decharge_1': 'decharge/fe_cal_out/free_ene.csv',
                'decharge_2': 'decharge_2/fe_cal_out/free_ene.csv',
                'decharge_3': 'decharge_3/fe_cal_out/free_ene.csv',
                'decharge_4': 'decharge_4/fe_cal_out/free_ene.csv',
                'decharge_5': 'decharge_5/fe_cal_out/free_ene.csv',
                'decharge_6': 'decharge_6/fe_cal_out/free_ene.csv',
                }
    csv_dict = {'5groups_plip': '/data/run01/scz1641/BACE1_bygroup/3rd_batch/lig_03_5group_plip_seq/openmm_run/complex/decharge/fe_cal_out/free_ene_group0.csv',
                '5groups_res_dis': '/data/run01/scz1641/BACE1_bygroup/3rd_batch/lig_03_5group_res_dist_based_seq/openmm_run/complex/decharge/fe_cal_out/free_ene_group0.csv',
                'com_decharge_1': '/HOME/scz1641/run/BACE1_bygroup/decharge_100wins_test/lig_03/openmm_run/complex/decharge/fe_cal_out/free_ene.csv',
                'com_decharge_2': '/HOME/scz1641/run/BACE1_bygroup/decharge_100wins_test/lig_03/openmm_run/complex/decharge_2/fe_cal_out/free_ene.csv',
                'com_decharge_3': '/HOME/scz1641/run/BACE1_bygroup/decharge_100wins_test/lig_03/openmm_run/complex/decharge_3/fe_cal_out/free_ene.csv',
                'com_decharge_4': '/HOME/scz1641/run/BACE1_bygroup/decharge_100wins_test/lig_03/openmm_run/complex/decharge_4/fe_cal_out/free_ene.csv',
                'com_decharge_5': '/HOME/scz1641/run/BACE1_bygroup/decharge_100wins_test/lig_03/openmm_run/complex/decharge_5/fe_cal_out/free_ene.csv',
                'com_decharge_6': '/HOME/scz1641/run/BACE1_bygroup/decharge_100wins_test/lig_03/openmm_run/complex/decharge_6/fe_cal_out/free_ene.csv',
                'lig_decharge_2': '/HOME/scz1641/run/BACE1_bygroup/decharge_100wins_test/lig_03/openmm_run/ligand/decharge_2/fe_cal_out/free_ene.csv',
                '9groups_plip': '/data/run01/scz1641/BACE1_bygroup/3rd_batch/lig_03_9group_plip_seq/openmm_run/complex/decharge/fe_cal_out/free_ene_group0.csv',
                '9groups_res_dis': '/data/run01/scz1641/BACE1_bygroup/3rd_batch/lig_03_9group_res_dist_based_seq/openmm_run/complex/decharge/fe_cal_out/free_ene_group0.csv',
                'normal_202_wins': '/data/run01/scz1641/BACE1_bygroup/3rd_batch/lig_03_normal_203_wins/openmm_run/complex/fe_cal_out/free_ene.csv',
                # 'decharge_3': 'decharge_3/fe_cal_out/free_ene.csv',
               }
    csv_dict = {'lig_decharge_1': '/HOME/scz1641/run/BACE1_bygroup/decharge_100wins_test/lig_03/openmm_run/ligand/decharge/fe_cal_out/free_ene.csv',
                'lig_decharge_2': '/HOME/scz1641/run/BACE1_bygroup/decharge_100wins_test/lig_03/openmm_run/ligand/decharge_2/fe_cal_out/free_ene.csv',
                'lig_decharge_3': '/HOME/scz1641/run/BACE1_bygroup/decharge_100wins_test/lig_03/openmm_run/ligand/decharge_3/fe_cal_out/free_ene.csv',
                'lig_5groups_plip': '/data/run01/scz1641/BACE1_bygroup/3rd_batch/lig_03_5group_plip_seq/openmm_run/ligand/decharge/fe_cal_out/free_ene_group0.csv',
                'lig_5groups_res_dis': '/data/run01/scz1641/BACE1_bygroup/3rd_batch/lig_03_5group_res_dist_based_seq/openmm_run/ligand/decharge/fe_cal_out/free_ene_group0.csv',
                'lig_9groups_plip': '/data/run01/scz1641/BACE1_bygroup/3rd_batch/lig_03_9group_plip_seq/openmm_run/ligand/decharge/fe_cal_out/free_ene_group0.csv',
                'lig_9groups_res_dis': '/data/run01/scz1641/BACE1_bygroup/3rd_batch/lig_03_9group_res_dist_based_seq/openmm_run/ligand/decharge/fe_cal_out/free_ene_group0.csv',
                'Normal_202_wins': '/data/run01/scz1641/BACE1_bygroup/3rd_batch/lig_03_normal_203_wins/openmm_run/ligand/fe_cal_out/free_ene.csv',
               }
    # ligs_lst = ['lig_02','lig_04','lig_05','lig_06', 'lig_07', 'lig_1a', 'lig_11', 'lig_36', 'lig_41', 'lig_45', 'lig_67', 'lig_13', 'lig_16', 'lig_32', 'lig_74', 'lig_81', 'lig_69']
    ligs_lst = ['lig_05']
    for i in ligs_lst:
        lig_name = i
        csv_dict = {f"{lig_name}_recharge_1": f'/HOME/scz1641/run/BACE1_bygroup/3rd_batch/{lig_name}/openmm_run/complex/decharge/fe_cal_out/free_ene_group0.csv',
                    f"{lig_name}_recharge_2": f'/HOME/scz1641/run/BACE1_bygroup/3rd_batch/{lig_name}/openmm_run/complex/decharge_2/fe_cal_out/free_ene_group0.csv',
                    f"{lig_name}_recharge_3": f'/HOME/scz1641/run/BACE1_bygroup/3rd_batch/{lig_name}/openmm_run/complex/decharge_3/fe_cal_out/free_ene_group0.csv',
                    f"{lig_name}_decharge_1": f'/HOME/scz1641/run/BACE1_bygroup/3rd_batch/{lig_name}/openmm_run/complex/decharge_an/fe_cal_out/free_ene_group0.csv',
                    f"{lig_name}_ligand_recharge_1": f'/HOME/scz1641/run/BACE1_bygroup/3rd_batch/{lig_name}/openmm_run/ligand/decharge/fe_cal_out/free_ene_group0.csv',
                    f"{lig_name}_recharge_4": f'/HOME/scz1641/run/BACE1_bygroup/3rd_batch/{lig_name}_8groups_manu/openmm_run/complex/decharge_2/fe_cal_out/free_ene_group0.csv',
                    f"{lig_name}_recharge_5": f'/HOME/scz1641/run/BACE1_bygroup/3rd_batch/{lig_name}_8groups_manu/openmm_run/complex/decharge_3/fe_cal_out/free_ene_group0.csv',
                    # f"{lig_name}_ligand_decharge_1": f'/HOME/scz1641/run/BACE1_bygroup/3rd_batch/{lig_name}/openmm_run/ligand/decharge_an/fe_cal_out/free_ene_group0.csv'
                    }
        ts = All_sys_data(csv_dict)
        ts.plot_dG_dlambda(f"{lig_name}_com_charges.png")
    
