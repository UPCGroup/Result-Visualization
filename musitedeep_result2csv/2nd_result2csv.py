import re
import pandas as pd

# 文件地址
result_txt_filename = '2nd_result2csv\\result_test_general_Y.txt'
result_csv_filename = '2nd_result2csv\\result_test_general_Y.csv'

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
            
            matchObj = re.match(r'\">sp\|(.*?)\|.*?\"\t(.*?)\t.*\"\t([.0-9]*)', str)
            if matchObj is None:
                print(str)
                print('###############ERROR################')
                    
            prots.append(matchObj.group(1))
            poss.append(matchObj.group(2))
            labels.append(matchObj.group(3))

        f.close()
    
    df = pd.DataFrame(data={'prot':prots, 'pos': poss, 'label': labels})
    df.to_csv(result_csv_filename, sep=',', index=False)
