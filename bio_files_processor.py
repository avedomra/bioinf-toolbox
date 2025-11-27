import os

def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str): 
    
    """
    Converts a multi-line FASTA format to one-line FASTA format.

    Arguments:
        input_fasta: path to the input FASTA file
        output_fasta: output file name (if not provided, will use the same name with '_oneline' suffix)
    """
    
    if not output_fasta:
        name, ext = os.path.splitext(input_fasta)
        output_fasta = f"{name}_oneline{ext if ext else '.fasta'}" 

    with open(input_fasta) as infile, open(output_fasta, "w") as outfile: 
        header = None
        seq = []

        for line in infile:
            line = line.strip()
            if not line:
                continue 

            if line.startswith(">"):
                if header:   
                    outfile.write(f"{header}\n{''.join(seq)}\n") 
                header = line
                seq = []
            else:
                seq.append(line) 

        if header:
            outfile.write(f"{header}\n{''.join(seq)}\n")

# это комментарии, которые я решила убрать из самой функции, но оставить внизу для себя, чтобы не сойти с ума :)
# сорри фо май рашн.... плиз игнор 

# если пользователь не передал аргумент "output_fasta", то он генерируется автоматически
# программа делит имя файла "input_fasta" и делит его на 'name' и 'extention'
# к 'name' добавляется '_oneline'
# к новому имени добавляется расширение инпут файла или '.fasta', если инпут без расширения

# цикл for - чтение инпут файла построчно 
# line.strip() - убирает пробельные символы
# if not line: continue - пропускает пустые строки

# если встретили заголовок (line.startswith(">"))
# if header - если заголовок заполнен, значит уже есть предыдущая запись
# добавляем предыдущий сиквенс в аутпут файл (outfile.write(f"{header}\n{''.join(seq)}\n"))
# теперь текущая строка становится новым заголовком, а seq - пустая строка
# если строка не начинается с ">" - добавляем ее в конец seq 

# цикл завершается на последней строке файла, но не записывает последний сиквенс
# потому что заканчивается на добавлении последней строки последнего сиквенса 
# поэтому нужно дозаписать после цикла 


def parse_blast_output(input_file: str, output_file: str):
    
    """
    Parses BLAST results and extracts names of the best matches for each query.

    Arguments:
        input_file (str): path to the input file
        output_file (str): path to the output file
    """
    
    proteins = []
    
    with open(input_file, 'r') as infile:
        content = infile.read()
    
    queries = content.split('Query #')

    for query in queries[1:]:
        lines = query.split('\n') 
        alignment_section = False 
        best_match_found = False 

        for line in lines:
            line = line.strip()

            if "Sequences producing significant alignments:" in line:
                alignment_section = True
                continue

            if alignment_section and line and not best_match_found:
                if line.startswith('-') or not line:
                    continue

            parts = [part for part in line.split('  ') if part.strip()]
            if parts:
                if '[' in line:
                    prot_name = line.split('[')[0].strip()
                else: 
                    prot_name = line.split('  ')[0].strip() if '  ' in line else line
                if prot_name.endswith('...'):
                    prot_name = prot_name[:-3]
                    proteins.append(prot_name)
                    best_match_found = True
    
    unique_proteins = sorted(set(proteins))

    with open(output_file, 'w') as outfile:
        for protein in unique_proteins:
            outfile.write(protein + '\n')


# задаем пустой список белков, которые будут записаны в аутпут файл 
# открываем инпут файл, переводим в строку 
# разделяем весь файл на запросы по ключевому выражению ("Query #") - это запросы по бласту для каждой последовательности

# запускаем большой цикл for, который пробегает по каждому запросу 
# сначала убираем заголовок всего txt файла [1:]
# делим запрос на строки, разделитель - перенос строки 
# задаем два состояния: 
# 1 - находимся ли мы в секции выравниваний 
# 2 - найдена ли там нужная последовательность (с лучшим мэтчем)
# в начале цикла - фолс для обоих состояний  

# запускаем маленький цикл for, который пробегает по каждой строке в запросе 
# доходим до заголовка "Sequences producing significant alignments" - теперь мы в секции выравниваний 
# меняем состояние на тру и пропускаем строку с названием
# задаем условие: if alignment_section and line and not best_match_found
# находимся в секции выравниваний, строка не пустая, нужная последовательность не найдена 
# теперь разбиваем строку на части по двойным пробелам 
# проверяем, записано ли название таксона в []
# если попадает, то указываем разделитель [, чтобы имя таксона не прилипло к имени белка 
# если нет, то продолжаем разделять по пробелам 
# если нет пробелов, то записываем строку целиком 
# убираем многоточие из названия белка 
# добавляем название белка в список белков 
# меняем состояние - нашли нужную последовательность 
# прерываем маленький цикл для строк внутри одного запроса, возвращаемся в большой цикл для разных запросов 

# делаем красивый список белков 
# оставляем только уникальные (set)
# сортируем в алфавитном порядке 
# записываем названия белков в аутпут файл с переносом строки в столбик 

# os.chdir()
# os.getcwd()
