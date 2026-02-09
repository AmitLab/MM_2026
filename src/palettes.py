pal_cytorisk = {'Standard Risk':(0.9427942677547513, 0.942825384792593, 0.9519953287278279),
                'Single Hit':(0.6285295118365009, 0.6285502565283954, 0.9679968858185519),
                '2+ Hits':(0.0, 0.0, 1.0),
                'NA':'#CCCCCC'}

pal_disease_old = {'AL': (0.12156862745098039, 0.4666666666666667, 0.7058823529411765),
 'EMD': (1.0, 0.4980392156862745, 0.054901960784313725),
 'Healthy': (0.17254901960784313, 0.6274509803921569, 0.17254901960784313),
 'MGUS': (0.8392156862745098, 0.15294117647058825, 0.1568627450980392),
 'MM_Unknown': (0.5803921568627451, 0.403921568627451, 0.7411764705882353),
 'NDAL': (0.5490196078431373, 0.33725490196078434, 0.29411764705882354),
 'NDMM': (0.8901960784313725, 0.4666666666666667, 0.7607843137254902),
 'PRMM': (0.4980392156862745, 0.4980392156862745, 0.4980392156862745),
 'RRMM': (0.7372549019607844, 0.7411764705882353, 0.13333333333333333),
 'SMM': (0.09019607843137255, 0.7450980392156863, 0.8117647058823529),
 'nan':'#CCCCCC',
 'non_naive_NDMM':'#a264af'}

pal_disease = {'Healthy':'#5bef8f',
            'MGUS':'#a0e3b7',
            'SMM':'#387472',
            'NDMM':'#c2c5e2',
            'non_naive_NDMM':'#887daf',
            'AL':'#3f436d',
            'RRMM':'#8c2e63',
            'MM_Unknown':'#CCCCCC',
            'NA':'#CCCCCC'}

pal_binary  = {'Yes':'#5f5f5f', 'No':'#FFFFFF', 'NA':'#CCCCCC'}

pal_binary_grey = {'Yes':'#181818', 'No':'#a9a9a9', 'NA':'#CCCCCC'}

pal_cnv_burden = {
    'High': '#E41A1C',   
    'Low': '#377EB8',     
    'NA':'#CCCCCC'
}

pal_cnv_call = {'ampl':'#EF3F3E',
                'neutral':'white',
                'del':'#4557A7',
                'NA':'#CCCCCC'}

pal_architype = {'5.0': (0.4980392156862745, 0.788235294117647, 0.4980392156862745),
 '3.0': '#15dec5',
 '2.0': '#2a8476',
 '7.0': '#91178b',
 #'nan': (0.2196078431372549, 0.4235294117647059, 0.6901960784313725),
 '6.0': '#5979fe',
 '4.0': '#2f4285',
 '8.0': (0.4, 0.4, 0.4),
 'NA': '#CCCCCC',
 '1.0': '#d7a547'}

pal_architype_renamed = {
 'MM2': '#15dec5',
 'MM1': '#2a8476',
 'MM5': '#91178b',
 'nan': '#CCCCCC',
  'None': '#CCCCCC',
  'NA': '#CCCCCC',
 'MM4': '#5979fe',
 'MM3': '#2f4285',
 'MM6':'#cdae4d'
 }

pal_cnv_call = {'ampl':'#e69892',
                'neutral':'white',
                'del':'#92b2e6',
                'NA':'#CCCCCC'}

pal_prolif_cnv_burden = {'Prolif_High':'#ce0c38', 'Prolif_Low':'#BB887C',
                         'Non-Prolif_High':'#89B2A5', 'Non-Prolif_Low':'#56DCCE'}

pal_prolif = {'True':'black', 'False':'grey', 'NA':'#CCCCCC'}

pal_cnv_call_grey = {'ampl': '#e69892', 'neutral': '#818181', 'del': '#92b2e6', 'NA': '#CCCCCC'}

pal_method = {'SPID':'#2A67C4',
              'MARS':'#2AC42F'}

pal_tissue = {'Blood':'#DB3838',
              'BM':'#3896DB',
              'EMD':'#38DB58',
              'NA':'#CCCCCC'
              }

pal_cell = {'Malignant': '#881d2a',
              'Healthy':'#2D8B57',
              'Healthy_Like':'#8EBC8F',
              'NA':'#CCCCCC'
              }

pal_cell_pb = {'Malignant':'#881d2a',
'Normal_Pb':'#d6ca72',
'Normal_PC':'#2D8B57',
'Interm':'#8EBC8F',
'NA':'#CCCCCC'}

pal_cell_renamed = {'Malignant': '#881d2a',
              'Normal PC':'#2D8B57',
              'Prolif. Normal':'#8EBC8F',
              'NA':'#CCCCCC'
              }

pal_trial = {'nan': '#CCCCCC',
 'KPT': (0.8509803921568627, 0.37254901960784315, 0.00784313725490196),
 'PPIA': (0.4588235294117647, 0.4392156862745098, 0.7019607843137254),
 'CART': (0.9058823529411765, 0.1607843137254902, 0.5411764705882353),
 'Transplantation': (0.9019607843137255,
  0.6705882352941176,
  0.00784313725490196),
 'NA': '#CCCCCC'}

pal_time = {'Pre': (0.5529411764705883, 0.8274509803921568, 0.7803921568627451),
 'Post': (1.0, 1.0, 0.7019607843137254),
 'C4': (0.7450980392156863, 0.7294117647058823, 0.8549019607843137),
 'C10': (0.984313725490196, 0.5019607843137255, 0.4470588235294118),
 'C1D22': (0.5019607843137255, 0.6941176470588235, 0.8274509803921568),
 'C3D1': (0.9921568627450981, 0.7058823529411765, 0.3843137254901961),
 'NA': '#CCCCCC'}

pal_iss = {'NA':'#CCCCCC',
           '1.0':'#D1EB42',
           '2.0':'#DBC738',
           '3.0':'#DB8238'}

pal_response = {'R': '#5885d5', 
                'NR': '#d44f79', 
                'NA':'#CCCCCC'}

pal_cells_tme = {
    'PC': '#1f77b4',
    'Mf': '#aec7e8',
    'T_Effector': '#ff7f0e',
    'T_Naive': '#ffbb78',
    'B': '#2ca02c',
    'Mast': '#98df8a',
    'T_Effector_GZMB': '#ff9896',
    'NK': '#9467bd',
    'B_Pro': '#c5b0d5',
    'pDC': '#8c564b',
    'Mo_CD16': '#c49c94',
    'Mo': '#e377c2',
    'Fibro': '#7f7f7f',
    'Neu_Pro': '#c7c7c7',
    'DC_IRF8': '#bcbd22',
    'DC': '#dbdb8d',
    'Mo_Pro': '#17becf',
    'Mega': '#9edae5',
    'Erythrocytes':'#b43b3b',
        'Malignant': '#881d2a', 
              'Healthy':'#2D8B57',
              'Healthy_Like':'#8EBC8F',
}

pal_cells_tme_with_pb = {
    'PC': '#1f77b4',
    'Mf': '#aec7e8',
    'T_Effector': '#ff7f0e',
    'T_Naive': '#ffbb78',
    'B': '#2ca02c',
    'Mast': '#98df8a',
    'T_Effector_GZMB': '#ff9896',
    'NK': '#9467bd',
    'B_Pro': '#c5b0d5',
    'pDC': '#8c564b',
    'Mo_CD16': '#c49c94',
    'Mo': '#e377c2',
    'Fibro': '#7f7f7f',
    'Neu_Pro': '#c7c7c7',
    'DC_IRF8': '#bcbd22',
    'DC': '#dbdb8d',
    'Mo_Pro': '#17becf',
    'Myeloid_Pro':'#17becf',
    'Mega': '#9edae5',
    'Erythrocytes':'#b43b3b',
        'Malignant': '#881d2a', 
              'Healthy':'#2D8B57',
              'Healthy_Like':'#8EBC8F',
              'Normal_Pb':'#d6ca72',
'Normal_PC':'#2D8B57',
'Interm':'#8EBC8F'
}

pal_cells_tme_wide = {
    'PC': '#1f77b4',
    'Mf': '#aec7e8',
    'T_Effector': '#ff7f0e',
    'T_Naive': '#ffbb78',
    'B': '#2ca02c',
    'Mast': '#98df8a',
    #'T_Effector_GZMB': '#ff9896',
    'NK': '#9467bd',
    'B_Pro': '#c5b0d5',
    'pDC': '#8c564b',
    #'Mo_CD16': '#c49c94',
    'Mo': '#e377c2',
    'Fibro': '#7f7f7f',
    #'Neu_Pro': '#c7c7c7',
    'Myeloid_Pro':'#17becf',
    'DC_IRF8': '#bcbd22',
    #'DC': '#dbdb8d',
    #'cDC1':'#dbdb8d',
    'cDC2':'#bcbd22',
    #'Mo_Pro': '#17becf',
    'Mega': '#9edae5',
    'Erythrocytes':'#b43b3b',
        'Malignant': '#881d2a', 
              'Healthy':'#2D8B57',
              'Healthy_Like':'#8EBC8F',
              'Normal_Pb':'#d6ca72',
'Normal_PC':'#2D8B57',
'Interm':'#8EBC8F'
}

pal_cells_tme_updated ={
'PC':'#7C7C7C',
'Malignant': '#7A1E1E',
'Normal_PC':'#55A279',
'Interm':'#8EBC8F',
'Normal_Pb':'#D6CA71',
###################
'T_Naive': '#93DCC9',
'T_Effector': '#318FA2',
'T_Effector_GZMB': '#01677F',
'NK':'#251188',
###################
'B':'#CC1F1A',
'B_Pro': '#F8BD47',
##################
'Fibro':'#E5EB40',
'Mega':'#4FE1F1',
'Mast':'#3798ED',
###################
'Mo_CD16':'#533F5E',
'Mf':'#ED59AD',
'DC':'#FAAFE3',
'Mo':'#827DB8',
'Mo_Pro': '#E2CBC1',
'Neu_Pro':'#D5D2E7',
'DC_IRF8': '#940393',
'pDC': '#81A6F2',
}

pal_cells_tme_updated_wide ={
'PC':'#7C7C7C',
'Malignant': '#7A1E1E',
'Normal_PC':'#55A279',
'Healthy':'#55A279',
'Interm':'#8EBC8F',
'Healthy_Like':'#8EBC8F',
'Normal_Pb':'#D6CA71',
'Plasmablastic':'#D6CA71',
###################
'T_Naive': '#93DCC9',
'T_Effector': '#318FA2',
#'T_Effector_GZMB': '#01677F',
'NK':'#251188',
###################
'B':'#CC1F1A',
'B_Pro': '#F8BD47',
##################
'Fibro':'#E5EB40',
'Mega':'#4FE1F1',
'Mast':'#3798ED',
###################
'Myeloid_Pro':'#E2CBC1',
'Mo_CD16':'#533F5E',
'Mf':'#ED59AD',
'DC':'#FAAFE3',
'Mo':'#827DB8',
'Mo_Pro': '#E2CBC1',
'Neu_Pro':'#D5D2E7',
'DC_IRF8': '#940393',
'pDC': '#81A6F2',
}

pal_architype_prolif = {'MM5_Non-Prolif':'#91178b',
'MM6_Non-Prolif':'#cdae4d',
                        'MM3_Non-Prolif':'#2f4285',
                        'MM4_Non-Prolif':'#5979fe',
       'Ambig_Prolif':'#797979',
       'MM1_Non-Prolif':'#2a8476',
       'MM5_Prolif':'#662d64',
       'MM6_Prolif':'#7a6936',
       'MM3_Prolif':'#424961',
       'MM4_Prolif':'#8c9ad3',
       'MM2_Non-Prolif':'#15dec5',
       'Ambig_Non-Prolif':'#b1b1b1',
       'Ambig':'#b1b1b1',
       'MM1_Prolif':'#1d4942',
       'MM2_Prolif':'#75aea7'
        }

pal_ms = {'CD1':'#1fbfd0',
'CD2':'#fdba13',
'HP':'#81c341',
'LB':'#8d574d',
'MF':'#d77ab1',
'MS':'#5852a3',
'PR':'#ed1d30',}
