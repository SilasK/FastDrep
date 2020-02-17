from snakemake.shell import shell

def run_bb(command,*args,log=None,**kwargs):
    for a in args:
        command+=' '+a

    for k in kwargs:
        if k=='input':
            input=kwargs[k]
            if len(input)>1 and not  isinstance(input,str):
                command+=f" in={','.join(kwargs[k])} "
            else:
                command+=f" in={kwargs[k]} "
        elif k=='mem':
            mem= int(float(kwargs[k])*0.85)
            command+=f" -Xmx{mem}g "
        else:
            if type(kwargs[k])==bool:
                kwargs[k]= 't' if kwargs[k] else 'f'
            command+=f" {k}={kwargs[k]} "

    if log is not None:
        if kwargs.get('append',False):
            command += f' 2>>{log}'
        else:
            command += f' 2>{log}'
    else:
        command += f' 2> /dev/null'


    #print(f"{command}")
    shell(command)


if __name__ == '__main__':

    if len(snakemake.log)>0:
        log=snakemake.log[0]
    else:
        log=None

    run_bb(log = log,
           **snakemake.input,
           **snakemake.params,
           **snakemake.output)
