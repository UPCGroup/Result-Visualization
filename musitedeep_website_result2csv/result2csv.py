import re
import pandas as pd

# 文件地址
result_txt_filename = 'musitedeep_result2csv\Prediction_results.txt'
result_csv_filename = 'musitedeep_result2csv\data_musitedeep.csv'

if __name__ == '__main__':

    # 处理后的数据信息
    prots = []
    poss = []
    rawseqs = []
    labels = []

    with open(result_txt_filename, 'r') as f:
        while True:
            str = f.readline()
            # 末尾行终止循环
            if str == '\n' or str=='':
                break
            # 跳过非数据行
            if len(list(str)) > 0 and list(str)[0] == '>' or str == 'ID	Position	Residue	PTMscores	Cutoff=0.5\n':
                continue
            
            matchObj = re.match(r'sp\|(.*?)\|.*?\t([0-9]*?)\t.*?Phosphoserine:(.*?)\t.*', str)
            if matchObj is None:
                matchObj = re.match(r'sp\|(.*?)\|.*?\t([0-9]*?)\t.*?Phosphothreonine:(.*?)\t.*', str)
            if matchObj is None:
                matchObj = re.match(r'sp\|(.*?)\|.*?\t([0-9]*?)\t.*?Phosphotyrosine:(.*?)\t.*', str)
            if matchObj is None:
                print(str)
                print('###############ERROR################')
                    
            prots.append(matchObj.group(1))
            poss.append(matchObj.group(2))
            labels.append(matchObj.group(3))

        f.close()
    
    df = pd.DataFrame(data={'prot':prots, 'pos': poss, 'label': labels})
    df.to_csv(result_csv_filename, sep=',', index=False)
