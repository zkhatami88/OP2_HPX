##########################################################################
#
# CUDA code generator
#
# This routine is called by op2 which parses the input files
#
# It produces a file xxx_kernel.cu for each kernel,
# plus a master kernel file
#
##########################################################################

import re
import datetime

def comm(line):
  global file_text, FORTRAN, CPP
  global depth
  prefix = ' '*depth
  if len(line) == 0:
    file_text +='\n'
  elif FORTRAN:
    file_text +='!  '+line+'\n'
  elif CPP:
    file_text +=prefix+'//'+line.rstrip()+'\n'

def rep(line,m):
  global dims, idxs, typs, indtyps, inddims
  if m < len(inddims):
    line = re.sub('INDDIM',str(inddims[m]),line)
    line = re.sub('INDTYP',str(indtyps[m]),line)

  line = re.sub('INDARG','ind_arg'+str(m),line)
  line = re.sub('DIM',str(dims[m]),line)
  line = re.sub('ARG','arg'+str(m),line)
  line = re.sub('TYP',typs[m],line)
  line = re.sub('IDX',str(int(idxs[m])),line)
  return line

def code(text):
  global file_text, FORTRAN, CPP, g_m
  global depth
  if text == '':
    prefix = ''
  else:
    prefix = ' '*depth
  file_text += prefix+rep(text,g_m).rstrip()+'\n'


def FOR(i,start,finish):
  global file_text, FORTRAN, CPP, g_m
  global depth
  if FORTRAN:
    code('do '+i+' = '+start+', '+finish+'-1')
  elif CPP:
    code('for ( int '+i+'='+start+'; '+i+'<'+finish+'; '+i+'++ ){')
  depth += 2

def FOR_INC(i,start,finish,inc):
  global file_text, FORTRAN, CPP, g_m
  global depth
  if FORTRAN:
    code('do '+i+' = '+start+', '+finish+'-1')
  elif CPP:
    code('for ( int '+i+'='+start+'; '+i+'<'+finish+'; '+i+' = '+i+'+='+inc+' ){')
  depth += 2

def ENDFOR():
  global file_text, FORTRAN, CPP, g_m
  global depth
  depth -= 2
  if FORTRAN:
    code('enddo')
  elif CPP:
    code('}')

def IF(line):
  global file_text, FORTRAN, CPP, g_m
  global depth
  if FORTRAN:
    code('if ('+line+') then')
  elif CPP:
    code('if ('+ line + ') {')
  depth += 2

def ENDIF():
  global file_text, FORTRAN, CPP, g_m
  global depth
  depth -= 2
  if FORTRAN:
    code('endif')
  elif CPP:
    code('}')

def op2_gen_cuda_simple_hyb(master, date, consts, kernels,sets):

  global dims, idxs, typs, indtyps, inddims
  global FORTRAN, CPP, g_m, file_text, depth

  OP_ID   = 1;  OP_GBL   = 2;  OP_MAP = 3;

  OP_READ = 1;  OP_WRITE = 2;  OP_RW  = 3;
  OP_INC  = 4;  OP_MAX   = 5;  OP_MIN = 6;

  accsstring = ['OP_READ','OP_WRITE','OP_RW','OP_INC','OP_MAX','OP_MIN' ]

  any_soa = 0
  for nk in range (0,len(kernels)):
    any_soa = any_soa or sum(kernels[nk]['soaflags'])

##########################################################################
#  create new kernel file
##########################################################################

  for nk in range (0,len(kernels)):

    name  = kernels[nk]['name']
    nargs = kernels[nk]['nargs']
    dims  = kernels[nk]['dims']
    maps  = kernels[nk]['maps']
    var   = kernels[nk]['var']
    typs  = kernels[nk]['typs']
    accs  = kernels[nk]['accs']
    idxs  = kernels[nk]['idxs']
    inds  = kernels[nk]['inds']
    soaflags = kernels[nk]['soaflags']

    ninds   = kernels[nk]['ninds']
    inddims = kernels[nk]['inddims']
    indaccs = kernels[nk]['indaccs']
    indtyps = kernels[nk]['indtyps']
    invinds = kernels[nk]['invinds']
    mapnames = kernels[nk]['mapnames']
    invmapinds = kernels[nk]['invmapinds']
    mapinds = kernels[nk]['mapinds']
    nmaps = 0

    if ninds > 0:
      nmaps = max(mapinds)+1

    vec =  [m for m in range(0,nargs) if int(idxs[m])<0 and maps[m] == OP_MAP]

    if len(vec) > 0:
      unique_args = [1];
      vec_counter = 1;
      vectorised = []
      new_dims = []
      new_maps = []
      new_vars = []
      new_typs = []
      new_accs = []
      new_idxs = []
      new_inds = []
      new_soaflags = []
      for m in range(0,nargs):
          if int(idxs[m])<0 and maps[m] == OP_MAP:
            if m > 0:
              unique_args = unique_args + [len(new_dims)+1]
            temp = [0]*(-1*int(idxs[m]))
            for i in range(0,-1*int(idxs[m])):
              temp[i] = var[m]
            new_vars = new_vars+temp
            for i in range(0,-1*int(idxs[m])):
              temp[i] = typs[m]
            new_typs = new_typs+temp
            for i in range(0,-1*int(idxs[m])):
              temp[i] = dims[m]
            new_dims = new_dims+temp
            new_maps = new_maps+[maps[m]]*int(-1*int(idxs[m]))
            new_soaflags = new_soaflags+[0]*int(-1*int(idxs[m]))
            new_accs = new_accs+[accs[m]]*int(-1*int(idxs[m]))
            for i in range(0,-1*int(idxs[m])):
              new_idxs = new_idxs+[i]
            new_inds = new_inds+[inds[m]]*int(-1*int(idxs[m]))
            vectorised = vectorised + [vec_counter]*int(-1*int(idxs[m]))
            vec_counter = vec_counter + 1;
          else:
            if m > 0:
              unique_args = unique_args + [len(new_dims)+1]
            new_dims = new_dims+[dims[m]]
            new_maps = new_maps+[maps[m]]
            new_accs = new_accs+[int(accs[m])]
            new_soaflags = new_soaflags+[soaflags[m]]
            new_idxs = new_idxs+[int(idxs[m])]
            new_inds = new_inds+[inds[m]]
            new_vars = new_vars+[var[m]]
            new_typs = new_typs+[typs[m]]
            vectorised = vectorised+[0]
      dims = new_dims
      maps = new_maps
      accs = new_accs
      idxs = new_idxs
      inds = new_inds
      var = new_vars
      typs = new_typs
      soaflags = new_soaflags;
      nargs = len(vectorised);

      for i in range(1,ninds+1):
        for index in range(0,len(inds)+1):
          if inds[index] == i:
            invinds[i-1] = index
            break
    else:
      vectorised = [0]*nargs
      unique_args = range(1,nargs+1)

    cumulative_indirect_index = [-1]*nargs;
    j = 0;
    for i in range (0,nargs):
      if maps[i] == OP_MAP:
        cumulative_indirect_index[i] = j
        j = j + 1
#
# set two logicals
#
    j = 0
    for i in range(0,nargs):
      if maps[i] == OP_MAP and accs[i] == OP_INC:
        j = i
    ind_inc = j > 0

    j = 0
    for i in range(0,nargs):
      if maps[i] == OP_GBL and accs[i] <> OP_READ:
        j = i
    reduct = j > 0

##########################################################################
#  start with CUDA kernel function
##########################################################################

    FORTRAN = 0;
    CPP     = 1;
    g_m = 0;
    file_text = ''
    depth = 0

    comm('user function')

    code('__device__')
    f = open(name+'.h', 'r')
    kernel_text = f.read()
    kernel_text = kernel_text.replace(name,name+'_gpu')
    file_text += kernel_text
    f.close()

    comm('')
    comm(' CUDA kernel function')

    if FORTRAN:
      code('subroutine op_cuda_'+name+'(')
    elif CPP:
      code('__global__ void op_cuda_'+name+'(')

    depth = 2

    for g_m in range(0,ninds):
      if (indaccs[g_m]==OP_READ):
        code('const INDTYP *__restrict INDARG,')
      else:
        code('INDTYP *__restrict INDARG,')

    if nmaps > 0:
      k = []
      for g_m in range(0,nargs):
        if maps[g_m] == OP_MAP and (not mapnames[g_m] in k):
          k = k + [mapnames[g_m]]
          code('const int *__restrict opDat'+str(invinds[inds[g_m]-1])+'Map, ')
    for g_m in range(0,nargs):
      if maps[g_m] == OP_ID:
        if accs[g_m] == OP_READ:
          code('const TYP *__restrict ARG,')
        else:
          code('TYP *ARG,')
    #for g_m in range(0,nargs):
      elif maps[g_m] == OP_GBL:
        if accs[g_m] == OP_INC or accs[g_m] == OP_MIN or accs[g_m] == OP_MAX:
          code('TYP *ARG,')
        elif accs[g_m] == OP_READ and dims[g_m].isdigit() and int(dims[g_m])==1:
          code('const TYP *ARG,')

    if ninds>0:
      if CPP:
        code('int    block_offset, ')
        code('int   *blkmap,       ')
        code('int   *offset,       ')
        code('int   *nelems,       ')
        code('int   *ncolors,      ')
        code('int   *colors,       ')
        code('int   nblocks,       ')
        code('int   set_size) {    ')
    else:
      code('int   set_size ) {')
      code('')

    for g_m in range(0,nargs):
      if maps[g_m]==OP_GBL and accs[g_m]<>OP_READ:
        code('TYP ARG_l[DIM];')
        if accs[g_m] == OP_INC:
          FOR('d','0','DIM')
          code('ARG_l[d]=ZERO_TYP;')
          ENDFOR()
        else:
          FOR('d','0','DIM')
          code('ARG_l[d]=ARG[d+blockIdx.x*DIM];')
          ENDFOR()
      elif maps[g_m]==OP_MAP and accs[g_m]==OP_INC:
        code('TYP ARG_l[DIM];')

    for m in range (1,ninds+1):
      g_m = m -1
      v = [int(inds[i]==m) for i in range(len(inds))]
      v_i = [vectorised[i] for i in range(len(inds)) if inds[i] == m]
      if sum(v)>1 and sum(v_i)>0: #check this sum(v_i)
        if indaccs[m-1] == OP_INC:
          ind = int(max([idxs[i] for i in range(len(inds)) if inds[i]==m])) + 1
          code('INDTYP *ARG_vec['+str(ind)+'] = {'); depth += 2;
          for n in range(0,nargs):
            if inds[n] == m:
              g_m = n
              code('ARG_l,')
          depth -= 2
          code('};')
        else:
          ind = int(max([idxs[i] for i in range(len(inds)) if inds[i]==m])) + 1
          code('INDTYP *ARG_vec['+str(ind)+'];')
#
# lengthy code for general case with indirection
#
    if ninds>0:
      code('')
      if ind_inc:
        code('__shared__ int    nelems2, ncolor;')

      code('__shared__ int    nelem, offset_b;')
      code('')
      code('extern __shared__ char shared[];')
      code('')
      IF('blockIdx.x+blockIdx.y*gridDim.x >= nblocks')
      code('return;')
      ENDIF()
      IF('threadIdx.x==0')
      code('')
      comm('get sizes and shift pointers and direct-mapped data')
      code('')
      code('int blockId = blkmap[blockIdx.x + blockIdx.y*gridDim.x  + block_offset];')
      code('')
      code('nelem    = nelems[blockId];')
      code('offset_b = offset[blockId];')
      code('')

      if ind_inc:
        code('nelems2  = blockDim.x*(1+(nelem-1)/blockDim.x);')
        code('ncolor   = ncolors[blockId];')
        code('')

      ENDIF()
      code('__syncthreads(); // make sure all of above completed')

      if ind_inc:
        FOR_INC('n','threadIdx.x','nelems2','blockDim.x')
        code('int col2 = -1;')
        k = []
        for g_m in range(0,nargs):
          if maps[g_m] == OP_MAP and (not mapinds[g_m] in k):
            k = k + [mapinds[g_m]]
            code('int map'+str(mapinds[g_m])+'idx;')
        IF('n<nelem')
        comm('initialise local variables')

        for g_m in range(0,nargs):
          if maps[g_m]==OP_MAP and accs[g_m]==OP_INC:
            FOR('d','0','DIM')
            code('ARG_l[d] = ZERO_TYP;')
            ENDFOR()
      else:
        FOR_INC('n','threadIdx.x','nelem','blockDim.x')
        k = []
        for g_m in range(0,nargs):
          if maps[g_m] == OP_MAP and (not mapinds[g_m] in k):
            k = k + [mapinds[g_m]]
            code('int map'+str(mapinds[g_m])+'idx;')


      k = []
      for g_m in range(0,nargs):
        if maps[g_m] == OP_MAP and (not mapinds[g_m] in k):
          k = k + [mapinds[g_m]]
          code('map'+str(mapinds[g_m])+'idx = opDat'+str(invmapinds[inds[g_m]-1])+'Map[n + offset_b + set_size * '+str(int(idxs[g_m]))+'];')


#
# simple alternative when no indirection
#
    else:
      code('')
      comm('process set elements')
      FOR_INC('n','threadIdx.x+blockIdx.x*blockDim.x','set_size','blockDim.x*gridDim.x')

#
# kernel call
#
    code('')
    comm('user-supplied kernel call')
    line = name+'_gpu('
    prefix = ' '*len(name)
    a = 0 #only apply indentation if its not the 0th argument
    indent =''
    for m in range (0, nargs):
      if a > 0:
        indent = '     '+' '*len(name)

      if maps[m] == OP_GBL:
        if accs[m] == OP_READ:
          line += rep(indent+'ARG,\n',m)
        else:
          line += rep(indent+'ARG_l,\n',m);
        a =a+1
      elif maps[m]==OP_MAP and  accs[m]==OP_INC:
        line += rep(indent+'ARG_l,\n',m)
        a =a+1
      elif maps[m]==OP_MAP:
        if soaflags[m]:
          line += rep(indent+'ind_arg'+str(inds[m]-1)+'+map'+str(mapinds[m])+'idx,'+'\n',m)
        else:
          line += rep(indent+'ind_arg'+str(inds[m]-1)+'+map'+str(mapinds[m])+'idx*DIM,'+'\n',m)
        a =a+1
      elif maps[m]==OP_ID:
        if ninds>0:
          if soaflags[m]:
            line += rep(indent+'ARG+(n+offset_b),\n',m)
          else:
            line += rep(indent+'ARG+(n+offset_b)*DIM,\n',m)
          a =a+1
        else:
          if soaflags[m]:
            line += rep(indent+'ARG+n,\n',m)
          else:
            line += rep(indent+'ARG+n*DIM,\n',m)
          a =a+1
      else:
        print 'internal error 1 '

    code(line[0:-2]+');') #remove final ',' and \n

#
# updating for indirect kernels ...
#
    if ninds>0:
      if ind_inc:
        code('col2 = colors[n+offset_b];')
        ENDIF()
        code('')
        comm('store local variables')
        code('')
        FOR('col','0','ncolor')
        IF('col2==col')

        for g_m in range(0,nargs):
          if maps[g_m] == OP_MAP and accs[g_m] == OP_INC:
            for d in range(0,int(dims[g_m])):
              code('ARG_l['+str(d)+'] += ind_arg'+str(inds[g_m]-1)+'['+str(d)+'+map'+str(mapinds[g_m])+'idx*DIM];')
        for g_m in range(0,nargs):
          if maps[g_m] == OP_MAP and accs[g_m] == OP_INC:
            for d in range(0,int(dims[g_m])):
              code('ind_arg'+str(inds[g_m]-1)+'['+str(d)+'+map'+str(mapinds[g_m])+'idx*DIM] = ARG_l['+str(d)+'];')

        ENDFOR()
        code('__syncthreads();')
        ENDFOR()
    ENDFOR()

#
# global reduction
#
    if reduct:
       code('')
       comm('global reductions')
       code('')
       for m in range (0,nargs):
         g_m = m
         if maps[m]==OP_GBL and accs[m]<>OP_READ:
           FOR('d','0','DIM')
           if accs[m]==OP_INC:
             code('op_reduction<OP_INC>(&ARG[d+blockIdx.x*DIM],ARG_l[d]);')
           elif accs[m]==OP_MIN:
             code('op_reduction<OP_MIN>(&ARG[d+blockIdx.x*DIM],ARG_l[d]);')
           elif accs[m]==OP_MAX:
             code('op_reduction<OP_MAX>(&ARG[d+blockIdx.x*DIM],ARG_l[d]);')
           else:
             print 'internal error: invalid reduction option'
             sys.exit(2);
           ENDFOR()
    depth -= 2
    code('}')
    code('')

##########################################################################
# then C++ stub function
##########################################################################

    code('')
    comm('GPU host stub function')
    code('void op_par_loop_'+name+'_gpu(char const *name, op_set set,')
    depth += 2

    for m in unique_args:
      g_m = m - 1
      if m == unique_args[len(unique_args)-1]:
        code('op_arg ARG){')
        code('')
      else:
        code('op_arg ARG,')

    code('using hpx::lcos::local::dataflow;')
    code('using hpx::util::unwrapped;')
    code('return dataflow(unwrapped([&name](op_set set,')

    for m in unique_args:
      g_m = m - 1
      if m == unique_args[len(unique_args)-1]:
        code('op_arg ARG){');
        code('')
      else:
        code('op_arg ARG,')

    for g_m in range (0,nargs):
      if maps[g_m]==OP_GBL:
        code('TYP*ARGh = (TYP *)ARG.data;')

    code('hpx::parallel::dynamic_chunk_size dcs(500);')

    code('int nargs = '+str(nargs)+';')
    code('op_arg args['+str(nargs)+'];')
    code('')

    #print vectorised

    for g_m in range (0,nargs):
      u = [i for i in range(0,len(unique_args)) if unique_args[i]-1 == g_m]
      if len(u) > 0 and vectorised[g_m] > 0:
        code('ARG.idx = 0;')
        code('args['+str(g_m)+'] = ARG;')

        v = [int(vectorised[i] == vectorised[g_m]) for i in range(0,len(vectorised))]
        first = [i for i in range(0,len(v)) if v[i] == 1]
        first = first[0]
        FOR('v','1',str(sum(v)))
        code('args['+str(g_m)+' + v] = op_arg_dat(arg'+str(first)+'.dat, v, arg'+\
        str(first)+'.map, DIM, "TYP", '+accsstring[accs[g_m]-1]+');')
        ENDFOR()
        code('')

      elif vectorised[g_m]>0:
        pass
      else:
        code('args['+str(g_m)+'] = ARG;')

#
# start timing
#
    code('')
    comm(' initialise timers')
    code('double cpu_t1, cpu_t2, wall_t1, wall_t2;')
    code('op_timing_realloc('+str(nk)+');')
    code('op_timers_core(&cpu_t1, &wall_t1);')
    code('OP_kernels[' +str(nk)+ '].name      = name;')
    code('OP_kernels[' +str(nk)+ '].count    += 1;')
    code('if (OP_kernels[' +str(nk)+ '].count==1) op_register_strides();')

    code('')

#
#   indirect bits
#
    if ninds>0:
      code('')
      code('int    ninds   = '+str(ninds)+';')
      line = 'int    inds['+str(nargs)+'] = {'
      for m in range(0,nargs):
        line += str(inds[m]-1)+','
      code(line[:-1]+'};')
      code('')

      IF('OP_diags>2')
      code('printf(" kernel routine with indirection: '+name+'\\n");')
      ENDIF()

      code('')
      comm('get plan')
      code('#ifdef OP_PART_SIZE_'+ str(nk))
      code('  int part_size = OP_PART_SIZE_'+str(nk)+';')
      code('#else')
      code('  int part_size = OP_part_size;')
      code('#endif')
      code('')
      code('int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);')

#
# direct bit
#
    else:
      code('')
      IF('OP_diags>2')
      code('printf(" kernel routine w/o indirection:  '+ name + '");')
      ENDIF()
      code('')
      code('op_mpi_halo_exchanges_cuda(set, nargs, args);')

    IF('set->size > 0')
    code('')
#    if any_soa:
#      code('int op2_stride_internal = set->size + set->exec_size + set->nonexec_size;')
#      #code('op_decl_const_char(1, "int", sizeof(int), (char *)&op2_stride, "op2_stride");')
#      code('cutilSafeCall(cudaMemcpyToSymbol(op2_stride , &op2_stride_internal, sizeof(int)));');
#      code('')

#
# kernel call for indirect version
#
    if ninds>0:
      code('op_plan *Plan = op_plan_get(name,set,part_size,nargs,args,ninds,inds);')
      code('')


#
# transfer constants
#
    g = [i for i in range(0,nargs) if maps[i] == OP_GBL and accs[i] == OP_READ]
    if len(g)>0:
      comm('transfer constants to GPU')
      code('int consts_bytes = 0;')
      for m in range(0,nargs):
        g_m = m
        if maps[m]==OP_GBL and accs[m]==OP_READ:
          code('consts_bytes += ROUND_UP(DIM*sizeof(TYP));')

      code('reallocConstArrays(consts_bytes);')
      code('consts_bytes = 0;')

      for m in range(0,nargs):
        if maps[m]==OP_GBL and accs[m]==OP_READ:
          g_m = m
          code('ARG.data   = OP_consts_h + consts_bytes;')
          code('ARG.data_d = OP_consts_d + consts_bytes;')
          FOR('d','0','DIM')
          code('((TYP *)ARG.data)[d] = ARGh[d];')
          ENDFOR()
          code('consts_bytes += ROUND_UP(DIM*sizeof(TYP));')
      code('mvConstArraysToDevice(consts_bytes);')
      code('')


#
# transfer global reduction initial data
#

    if ninds == 0:
      comm('set CUDA execution parameters')
      code('#ifdef OP_BLOCK_SIZE_'+str(nk))
      code('  int nthread = OP_BLOCK_SIZE_'+str(nk)+';')
      code('#else')
      code('  int nthread = OP_block_size;')
      comm('  int nthread = 128;')
      code('#endif')
      code('')
      code('int nblocks = 200;')
      code('')

    if reduct:
      comm('transfer global reduction data to GPU')
      if ninds>0:
        code('int maxblocks = 0;')
        FOR('col','0','Plan->ncolors')
        code('maxblocks = MAX(maxblocks,Plan->ncolblk[col]);')
        ENDFOR()
      else:
        code('int maxblocks = nblocks;')

      code('int reduct_bytes = 0;')
      code('int reduct_size  = 0;')

      for g_m in range(0,nargs):
        if maps[g_m]==OP_GBL and accs[g_m]<>OP_READ:
          code('reduct_bytes += ROUND_UP(maxblocks*DIM*sizeof(TYP));')
          code('reduct_size   = MAX(reduct_size,sizeof(TYP));')

      code('reallocReductArrays(reduct_bytes);')
      code('reduct_bytes = 0;')

      for g_m in range(0,nargs):
        if maps[g_m]==OP_GBL and accs[g_m]<>OP_READ:
          code('ARG.data   = OP_reduct_h + reduct_bytes;')
          code('ARG.data_d = OP_reduct_d + reduct_bytes;')
          FOR('b','0','maxblocks')
          FOR('d','0','DIM')
          if accs[g_m]==OP_INC:
            code('((TYP *)ARG.data)[d+b*DIM] = ZERO_TYP;')
          else:
            code('((TYP *)ARG.data)[d+b*DIM] = ARGh[d];')
          ENDFOR()
          ENDFOR()
          code('reduct_bytes += ROUND_UP(maxblocks*DIM*sizeof(TYP));')
      code('mvReductArraysToDevice(reduct_bytes);')
      code('')

#
# kernel call for indirect version
#
    if ninds>0:
      comm('execute plan')
      code('')
      code('int block_offset = 0;')
      FOR('col','0','Plan->ncolors')
      IF('col==Plan->ncolors_core')
      code('op_mpi_wait_all_cuda(nargs, args);')
      ENDIF()
      code('#ifdef OP_BLOCK_SIZE_'+str(nk))
      code('int nthread = OP_BLOCK_SIZE_'+str(nk)+';')
      code('#else')
      code('int nthread = OP_block_size;')
      code('#endif')
      code('')
      code('dim3 nblocks = dim3(Plan->ncolblk[col] >= (1<<16) ? 65535 : Plan->ncolblk[col],')
      code('Plan->ncolblk[col] >= (1<<16) ? (Plan->ncolblk[col]-1)/65535+1: 1, 1);')
      IF('Plan->ncolblk[col] > 0')

      if reduct:
        code('int nshared = reduct_size*nthread;')
        code('op_cuda_'+name+'<<<nblocks,nthread,nshared>>>(')
      else:
        code('op_cuda_'+name+'<<<nblocks,nthread>>>(')

      for m in range(1,ninds+1):
        g_m = invinds[m-1]
        code('(TYP *)ARG.data_d,')
      if nmaps > 0:
        k = []
        for g_m in range(0,nargs):
          if maps[g_m] == OP_MAP and (not mapnames[g_m] in k):
            k = k + [mapnames[g_m]]
            code('arg'+str(invinds[inds[g_m]-1])+'.map_data_d, ')
      for g_m in range(0,nargs):
        if inds[g_m]==0:
          code('(TYP*)ARG.data_d,')


      code('block_offset,')
      code('Plan->blkmap,')
      code('Plan->offset,')
      code('Plan->nelems,')
      code('Plan->nthrcol,')
      code('Plan->thrcol,')
      code('Plan->ncolblk[col],')
      code('set->size+set->exec_size);')
      code('')
      if reduct:
        comm('transfer global reduction data back to CPU')
        IF('col == Plan->ncolors_owned-1')
        code('mvReductArraysToHost(reduct_bytes);')
        ENDIF()

      ENDFOR()
      code('block_offset += Plan->ncolblk[col];')
      ENDIF()
#
# kernel call for direct version
#
    else:
      if reduct:
        code('int nshared = reduct_size*nthread;')
        code('op_cuda_'+name+'<<<nblocks,nthread,nshared>>>(')
      else:
        code('op_cuda_'+name+'<<<nblocks,nthread>>>(')

      indent = '  '#*(len(name)+42)
      for g_m in range(0,nargs):
        if g_m > 0:
          code(indent+'(TYP *) ARG.data_d,')
        else:
          code(indent+'(TYP *) ARG.data_d,')

      code(indent+'set->size );')

    if ninds>0:
      code('OP_kernels['+str(nk)+'].transfer  += Plan->transfer;')
      code('OP_kernels['+str(nk)+'].transfer2 += Plan->transfer2;')


#
# transfer global reduction initial data
#
    if reduct:
      if ninds == 0:
        comm('transfer global reduction data back to CPU')
        code('mvReductArraysToHost(reduct_bytes);')

      for m in range(0,nargs):
        g_m = m
        if maps[m]==OP_GBL and accs[m]<>OP_READ:
          FOR('b','0','maxblocks')
          FOR('d','0','DIM')
          if accs[m]==OP_INC:
            code('ARGh[d] = ARGh[d] + ((TYP *)ARG.data)[d+b*DIM];')
          elif accs[m]==OP_MIN:
            code('ARGh[d] = MIN(ARGh[d],((TYP *)ARG.data)[d+b*DIM]);')
          elif accs[m]==OP_MAX:
            code('ARGh[d] = MAX(ARGh[d],((TYP *)ARG.data)[d+b*DIM]);')
          ENDFOR()
          ENDFOR()

          code('ARG.data = (char *)ARGh;')
          code('op_mpi_reduce(&ARG,ARGh);')

    ENDIF()
    code('op_mpi_set_dirtybit_cuda(nargs, args);')

#
# update kernel record
#

    code('cutilSafeCall(cudaDeviceSynchronize());')
    comm('update kernel record')
    code('op_timers_core(&cpu_t2, &wall_t2);')
    code('OP_kernels[' +str(nk)+ '].time     += wall_t2 - wall_t1;')

    if ninds == 0:
      line = 'OP_kernels['+str(nk)+'].transfer += (float)set->size *'

      for g_m in range (0,nargs):
        if maps[g_m]<>OP_GBL:
          if accs[g_m]==OP_READ or accs[g_m]==OP_WRITE:
            code(line+' ARG.size;')
          else:
            code(line+' ARG.size * 2.0f;')

    depth = depth - 2
    code('}),set,')

    for m in unique_args:
      g_m = m - 1
      if m == unique_args[len(unique_args)-1]:
        code('ARG');
        code('')
      else:
        code('ARG,')

    code(');')

    code('}')


    code('')
    code('void op_par_loop_'+name+'_cpu(char const *name, op_set set,')
    depth += 2

    for m in unique_args:
      g_m = m - 1
      if m == unique_args[len(unique_args)-1]:
        code('op_arg ARG);')
        code('')
      else:
        code('op_arg ARG,')
    depth -=2


    code('')
    comm('GPU host stub function')
    code('#if OP_HYBRID_GPU')
#    code('void op_par_loop_'+name+'(char const *name, op_set set,')
    code('hpx::future<void> op_par_loop_'+name+'(char const *name, hpx::future<op_set> set,')
    depth += 2

    for m in unique_args:
      g_m = m - 1
      if m == unique_args[len(unique_args)-1]:
        code('hpx::future<op_arg> ARG){')
        code('')
      else:
        code('hpx::future<op_arg> ARG,')

    IF('OP_hybrid_gpu')
    code('op_par_loop_'+name+'_gpu(name, set,')
    depth += 2
    for m in unique_args:
      g_m = m - 1
      if m == unique_args[len(unique_args)-1]:
        code('ARG);')
        code('')
      else:
        code('ARG,')
    depth -=2
    code('}else{')
    code('op_par_loop_'+name+'_cpu(name, set,')
    depth += 2
    for m in unique_args:
      g_m = m - 1
      if m == unique_args[len(unique_args)-1]:
        code('ARG);')
        code('')
      else:
        code('ARG,')
    depth -=2
    ENDIF()
    depth-=2
    code('}')
    code('#else')
    code('void op_par_loop_'+name+'(char const *name, op_set set,')
    depth += 2

    for m in unique_args:
      g_m = m - 1
      if m == unique_args[len(unique_args)-1]:
        code('op_arg ARG){')
        code('')
      else:
        code('op_arg ARG,')


    code('op_par_loop_'+name+'_gpu(name, set,')
    depth += 2
    for m in unique_args:
      g_m = m - 1
      if m == unique_args[len(unique_args)-1]:
        code('ARG);')
        code('')
      else:
        code('ARG,')
    depth-=2
    code('}')
    depth-=2
    code('#endif //OP_HYBRID_GPU')

##########################################################################
#  output individual kernel file
##########################################################################
    fid = open(name+'_kernel.cu','w')
    date = datetime.datetime.now()
    fid.write('//\n// auto-generated by op2.py\n//\n\n')
    fid.write(file_text)
    fid.close()

# end of main kernel call loop


##########################################################################
#  output one master kernel file
##########################################################################

  file_text = ''
  comm('header')
  code('#include "op_lib_cpp.h"')
  code('#include "op_cuda_rt_support.h"')
  code('#include "op_cuda_reduction.h"')
  code('')
  comm('global constants')
  code('#ifndef MAX_CONST_SIZE')
  code('#define MAX_CONST_SIZE 128')
  code('#endif')
  code('')

  code('#define STRIDE(x,y) x*y')
  for ns in range (0,len(sets)):
    code('__constant__ int '+sets[ns]['name'].replace('"','')+'_stride;')

  for nc in range (0,len(consts)):
    if consts[nc]['dim']==1:
      code('__constant__ '+consts[nc]['type'][1:-1]+' '+consts[nc]['name'].replace('"','')+';')
    else:
      if consts[nc]['dim'] > 0:
        num = str(consts[nc]['dim'])
      else:
        num = 'MAX_CONST_SIZE'

      code('__constant__ '+consts[nc]['type'][1:-1]+' '+consts[nc]['name'].replace('"','')+'['+num+'];')

  if any_soa:
    code('__constant__ int op2_stride;')
    code('')
    code('#define OP2_STRIDE(arr, idx) arr[op2_stride*(idx)]')

  code('')
  code('void op_register_strides() {')
  depth = depth + 2
  code('int size;')
  for ns in range (0,len(sets)):
    code('size = op_size_of_set("'+sets[ns]['name'].replace('"','')+'");')
    code('cutilSafeCall(cudaMemcpyToSymbol('+sets[ns]['name'].replace('"','')+'_stride, &size, sizeof(int)));')
  depth = depth - 2
  code('}')
  code('')
  code('void op_decl_const_char(int dim, char const *type,')
  code('int size, char *dat, char const *name){')
  depth = depth + 2
  IF('OP_hybrid_gpu')

  for nc in range(0,len(consts)):
    IF('!strcmp(name,"'+consts[nc]['name'].replace('"','')+'")')
    if consts[nc]['dim'] < 0:
      IF('!strcmp(name,"'+consts[nc]['name']+'") && size>MAX_CONST_SIZE) {')
      code('printf("error: MAX_CONST_SIZE not big enough\n"); exit(1);')
      ENDIF()
    code('cutilSafeCall(cudaMemcpyToSymbol('+consts[nc]['name'].replace('"','')+', dat, dim*size));')
    ENDIF()
    code('else ')

  code('{')
  depth = depth + 2
  code('printf("error: unknown const name\\n"); exit(1);')
  ENDIF()
  ENDIF()

  depth = depth - 2
  code('}')
  code('')
  comm('user kernel files')

  for nk in range(0,len(kernels)):
    file_text = file_text +\
    '#include "'+kernels[nk]['name']+'_kernel.cu"\n'

  master = master.split('.')[0]
  fid = open(master.split('.')[0]+'_kernels.cu','w')
  fid.write('//\n// auto-generated by op2.py\n//\n\n')
  fid.write(file_text)
  fid.close()

  fid= open(master.split('.')[0]+'_kernels.cpp','r')
  text = fid.read()
  text = text.replace('_kernel.cpp','_cpu_kernel.cpp')
  fid.close()
  fid= open(master.split('.')[0]+'_cpu_kernels.cpp','w')
  fid.write(text)
  fid.close()
  for nk in range (0,len(kernels)):
    name  = kernels[nk]['name']
    fid = open(name+'_kernel.cpp','r')
    text = fid.read()
    fid.close()
    fid = open(name+'_cpu_kernel.cpp','w')
    fid.write(text.replace('op_par_loop_'+name+'(', 'op_par_loop_'+name+'_cpu('))
    fid.close()
