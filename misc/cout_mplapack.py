from __future__ import absolute_import, division, print_function
from six.moves import range
from libtbx.utils import product
from libtbx import group_args
from libtbx import mutable
from libtbx import Auto
import os.path
from six.moves import zip
import re

def delete_brackets(s):
  l = [ '\([^()]+\)', '\[[^\[\]]+\]', '\[\]' , '\(\)' ]
  for _l in l:
    s = re.sub(_l, "", s)
  return delete_brackets(s) if sum([1 if re.search(_l, s) else 0 for _l in l]) > 0 else s
																 
fmt_comma_placeholder = chr(255)

def break_line_if_necessary(callback, line, max_len=80, min_len=70):
  def cb_finalize(line):
    callback(line.replace(fmt_comma_placeholder, ","))
  nc = len(line)
  if (nc <= max_len):
    cb_finalize(line)
    return
  for i_start in range(nc):
    if (line[i_start] != " "):
      break
  else:
    raise AssertionError
  lsw = line.startswith
  if (lsw("//", i_start)):
    cb_finalize(line)
    return
  potential_break_points = []
  ic = i_start
  while (ic < nc):
    c = line[ic]
    if ("\"'".find(c) >= 0):
      q = c
      ic_q = ic
      ic += 1
      while (ic < nc):
        prev_c = c
        c = line[ic]
        ic += 1
        if (c == q and prev_c != "\\"):
          break
      else:
        raise AssertionError
    elif (c == "("):
      ic += 1
      potential_break_points.append((0, ic))
    elif (lsw(", ", ic)):
      ic += 2
      potential_break_points.append((1, ic))
    elif (   lsw(" = ", ic)
          or lsw("), ", ic)):
      ic += 3
      potential_break_points.append((1, ic))
    elif (   lsw(" + ", ic)
          or lsw(" - ", ic)
          or lsw(" * ", ic)
          or lsw(" / ", ic)):
      ic += 3
      potential_break_points.append((0, ic))
    elif (   lsw(" && ", ic)
          or lsw(" || ", ic)):
      ic += 4
      potential_break_points.append((0, ic))
    elif (lsw("//", ic)):
      break
    else:
      ic += 1
  potential_break_points.append((0, nc))
  n = nc - i_start
  from libtbx.math_utils import iround, iceil
  l = max(min_len, iround(n / iceil(n / (max_len - i_start - 2))))
  b = 0
  f = 0
  def break_more_if_necessary(s):
    while (f+len(s) > max_len and s.startswith('"')):
      i = max_len-2-f
      j = s.rfind(fmt_comma_placeholder, 0, i)
      if (j > 4): # ad-hoc value
        i = j+1
      else:
        for j in range(i-1,-1,-1):
          if (s[j] != "\\"):
            if ((i - j ) % 2 == 0):
              i -= 1
            break
        else:
          raise AssertionError
      cb_finalize(" "*f + s[:i] + '"')
      s = '"' + s[i:]
    cb_finalize(" "*f + s)
  if (lsw("if (", i_start)): indent_width = 4
  else:                      indent_width = 2
  pprio = 0
  pp = 0
  for ip in range(len(potential_break_points)):
    prio,p = potential_break_points[ip]
    def following_point_is_better():
      for jp in range(ip,len(potential_break_points)):
        prio,p = potential_break_points[jp]
        if (prio == 1 and p-b+f <= max_len):
          return True
      return False
    if (    p-b+f > l
        and b != pp
        and (pprio == 1 or not following_point_is_better())):
      s = line[b:pp].rstrip()
      if (f == 0):
        cb_finalize(s)
        f = i_start + indent_width
      else:
        break_more_if_necessary(s=s)
      b = pp
    pprio = prio
    pp = p
  if (b < nc):
    break_more_if_necessary(s=line[b:])

def break_lines(cpp_text, prev_line=None):
  prev_line = [prev_line]
  result = []
  def callback(line):
    if (   prev_line[0] is None
        or line != prev_line[0]
        or not line.lstrip().startswith("//C")):
      result.append(line)
      prev_line[0] = line
  for line in "\n".join(cpp_text).splitlines():
    break_line_if_necessary(callback=callback, line=line)
  return result

class dynamic_parameter_props(object):

  __slots__ = ["name", "ctype", "default"]

  def __init__(O, name, ctype, default):
    O.name = name
    O.ctype = ctype
    O.default = default

def create_buffer_blocks(
      target_number_of_blocks,
      buffers,
      min_lines_per_block=100):
  if (target_number_of_blocks >= len(buffers)):
    return [[buffer] for buffer in buffers]
  numbers_of_lines = [len(lines) for lines in buffers]
  sum_lines = sum(numbers_of_lines)
  lines_per_block = max(
    min_lines_per_block,
    sum_lines / target_number_of_blocks)
  result = []
  block_sum_lines = 0
  j = 0
  for i,n in enumerate(numbers_of_lines):
    if (block_sum_lines + n > lines_per_block):
      if (i == j or block_sum_lines <= min_lines_per_block):
        result.append(buffers[j:i+1])
        block_sum_lines = 0
        j = i+1
      else:
        result.append(buffers[j:i])
        block_sum_lines = n
        j = i
    else:
      block_sum_lines += n
  if (j < len(buffers)):
    result.append(buffers[j:])
  assert sum([len(block) for block in result]) == len(buffers)
  return result

def show_traceback():
  import traceback
  print(traceback.format_exc(limit=None))

def strip_leading_zeros(string):
  for i in range(len(string)):
    if (string[i] != "0"):
      return string[i:]
  if (len(string) == 0):
    return ""
  return "0"

def escape_string_literal(s):
  return (s
    .replace("\\","\\\\")
    .replace('"','\\"')
    .replace("\t", "\\t")
    .replace("??", "\\?\\?"))

def convert_complex_literal(vmap, tok):
  assert len(tok.value) == 4
  cc = []
  for part in tok.value[1:3]:
    c = []
    sign_tok, val_tok = part
    if (sign_tok is not None):
      c.append(sign_tok.value)
    c.append(convert_token(vmap=vmap, leading=None, tok=val_tok))
    cc.append("".join(c))
  return "fem::cmplx(%s)" % ", ".join(cc)

def convert_token(vmap, leading, tok, had_str_concat=None):
  tv = tok.value
  if (tok.is_identifier()):
    return vmap.get(tv, tv)
  if (tok.is_op()):
    if (tv == ".not."):
      return "!"
    if (tv == ".and."):
      return " && "
    if (tv == ".or."):
      return " || "
    if (tv == ".eqv."):
      return " == "
    if (tv == ".neqv."):
      return " != "
    if (tv in ["+", "-"]):
      if (leading):
        return tv
      return " "+tv+" "
    if (tv == "*"):
      if (leading):
        return "star "
      return " "+tv+" "
    if (tv == "/"):
      return " "+tv+" "
    if (tv == "//"):
      if (had_str_concat is None): tok.raise_internal_error()
      had_str_concat.value = True
      return " + "
    if (tv == ":"):
      if (leading):
        return "1, "
      return ", "
    if (tv == ".eq." or tv == "=="):
      return " == "
    if (tv == ".ne." or tv == "/="):
      return " != "
    if (tv == ".lt." or tv == "<"):
      return " < "
    if (tv == ".le." or tv == "<="):
      return " <= "
    if (tv == ".gt." or tv == ">"):
      return " > "
    if (tv == ".ge." or tv == ">="):
      return " >= "
    tok.raise_not_supported()
  if (tok.is_string()):
    s = '"' + escape_string_literal(tok.value) + '"'
    if (had_str_concat is None or not had_str_concat.value):
      return s
    return "str_cref(%s)" % s
  if (tok.is_logical()):
    if (tv == ".false."):
      return "false"
    return "true"
  if (tok.is_integer()):
    return tv
  if (tok.is_hexadecimal()):
    return "0x"+tv
  if (tok.is_real()):
    return tv+"f"
  if (tok.is_double_precision()):
    return tv.replace("d", "e")
  if (tok.is_complex()):
    return convert_complex_literal(vmap=vmap, tok=tok)
  tok.raise_not_supported()

class major_types_cache(object):

  __slots__ = ["identifiers"]

  def __init__(O):
    O.identifiers = None

  def __contains__(O, value):
    if (O.identifiers is None):
      O.identifiers = set()
      import libtbx.load_env
      hpp = libtbx.env.under_dist(
        module_name="fable", path="fem/major_types.hpp", test=os.path.isfile)
      using_fem = "  using fem::"
      with open(hpp) as f:
        lines = f.read().splitlines()
      for line in lines:
        if (line.startswith(using_fem)):
          assert line.endswith(";")
          O.identifiers.add(line[len(using_fem):-1])
    return value in O.identifiers

major_types = major_types_cache()

cpp_keywords = set("""\
and and_eq asm auto bitand bitor bool break case catch char class compl const
const_cast continue default delete do double dynamic_cast else enum explicit
export extern false float for friend goto if inline int long mutable namespace
new not not_eq operator or or_eq private protected public register
reinterpret_cast return short signed sizeof static static_cast struct switch
template this throw true try typedef typeid typename union unsigned using
virtual void volatile wchar_t while xor xor_eq argv argc
""".split())

def prepend_identifier_if_necessary(identifier):
  if (identifier in major_types or identifier in cpp_keywords):
    return "identifier_" + identifier
  return identifier

def produce_comment_given_sl(callback, sl):
  if (sl.stmt_offs is None):
    t = sl.text[1:]
  elif (sl.index_of_exclamation_mark is not None):
    t = sl.stmt[sl.index_of_exclamation_mark+1:]
  else:
    t = None
  if (t is not None):
    callback("//C%s" % t.expandtabs().rstrip())

def produce_comment_given_ssl(callback, ssl):
  if (ssl is None): return
  for sl in ssl.source_line_cluster:
    produce_comment_given_sl(callback=callback, sl=sl)
  return

def produce_comments(callback, ssl_list):
  for ssl in ssl_list:
    produce_comment_given_ssl(callback=callback, ssl=ssl)

def flush_comments_if_non_trivial(callback, buffer):
  for line in buffer:
    if (line != "//C"):
      for line in buffer:
        callback(line)
      return

def produce_leading_comments(callback, fproc):
  buffer = []
  produce_comments(callback=buffer.append, ssl_list=fproc.leading_comments)
  produce_comment_given_ssl(
    callback=buffer.append, ssl=fproc.top_ssl)
  flush_comments_if_non_trivial(callback=callback, buffer=buffer)

def produce_trailing_comments(callback, fproc):
  buffer = []
  produce_comment_given_ssl(callback=buffer.append, ssl=fproc.end_ssl)
  produce_comments(
    callback=buffer.append, ssl_list=fproc.trailing_comments)
  flush_comments_if_non_trivial(callback=callback, buffer=buffer)

class comment_manager(object):

  __slots__ = ["sl_list", "index"]

  def __init__(O, fproc):
    O.sl_list = []
    for ssl in fproc.body_lines:
      if (ssl is not None):
        for sl in ssl.source_line_cluster:
          O.sl_list.append(sl)
    O.sl_list.sort(key=lambda source_line: source_line.global_line_index)
    O.index = 0

  def produce(O, callback):
    produce_comment_given_sl(callback=callback, sl=O.sl_list[O.index])
    O.index += 1

  def insert_before(O, executable_info, callback):
    i = executable_info.ssl.source_line_cluster[-1].global_line_index
    while (O.index != len(O.sl_list)):
      j = O.sl_list[O.index].global_line_index
      if (j > i):
        break
      O.produce(callback=callback)

  def flush_remaining(O, callback):
    while (O.index != len(O.sl_list)):
      O.produce(callback=callback)

class conv_hook_info(object):

  __slots__ = [
    "ignore_common_and_save",
    "needs_sve_dynamic_parameters",
    "variant_common_names",
    "needs_is_called_first_time",
    "needs_variant_bind",
    "data_init_after_variant_bind"]

  def __init__(O):
    O.ignore_common_and_save = False
    O.needs_sve_dynamic_parameters = False
    O.variant_common_names = None
    O.needs_is_called_first_time = None
    O.needs_variant_bind = None
    O.data_init_after_variant_bind = None

class global_conversion_info(object):

  __slots__ = [
    "topological_fprocs",
    "dynamic_parameters",
    "fortran_file_comments",
    "fem_do_safe",
    "arr_nd_size_max",
    "inline_all",
    "fprocs_by_name",
    "converted_commons_info",
    "separate_namespaces",
    "data_values_block_size",
    "data_specializations"]

  def __init__(O,
        topological_fprocs,
        dynamic_parameters,
        fortran_file_comments,
        fem_do_safe,
        arr_nd_size_max,
        inline_all,
        converted_commons_info,
        separate_namespaces,
        data_values_block_size,
        data_specializations):
    O.topological_fprocs = topological_fprocs
    O.dynamic_parameters = dynamic_parameters
    O.fortran_file_comments = fortran_file_comments
    O.fem_do_safe = fem_do_safe
    O.arr_nd_size_max = arr_nd_size_max
    O.inline_all = inline_all
    O.fprocs_by_name = topological_fprocs.all_fprocs.fprocs_by_name()
    O.converted_commons_info = converted_commons_info
    O.separate_namespaces = separate_namespaces
    O.data_values_block_size = data_values_block_size
    O.data_specializations = data_specializations

  def specialized(O, fproc):
    return conversion_info(global_conv_info=O, fproc=fproc)

class conversion_info(global_conversion_info):

  __slots__ = global_conversion_info.__slots__ + [
    "fproc",
    "comment_manager",
    "vmap"]

  def __init__(O,
        global_conv_info=None,
        fproc=None,
        vmap=None):
    val = None
    for slot in global_conversion_info.__slots__:
      if (global_conv_info is not None):
        val = getattr(global_conv_info, slot)
      setattr(O, slot, val)
    O.fproc = fproc
    if (O.fproc is None):
      O.comment_manager = None
    else:
      O.comment_manager = comment_manager(fproc=O.fproc)
    if (vmap is None):
      O.vmap = {}
    else:
      O.vmap = vmap

  def set_vmap_force_local(O, fdecl):
    identifier = fdecl.id_tok.value
    O.vmap[identifier] = prepend_identifier_if_necessary(identifier)

  def set_vmap_for_callable(O, identifier):
    if (O.separate_namespaces is not None):
      ns = O.separate_namespaces.get(identifier)
      if (ns is not None):
        O.vmap[identifier] \
          = ns + "::" + prepend_identifier_if_necessary(identifier)
        return True
    if (identifier in ["getargc", "iargc"]):
      O.vmap[identifier] = "cmn.%s" % identifier
      return True
    from fable import intrinsics
    if (identifier in intrinsics.extra_set_lower):
      O.vmap[identifier] = "fem::" + identifier
      return True
    return False

  def set_vmap_from_fdecl(O, fdecl):
    identifier = fdecl.id_tok.value
    if (fdecl.is_common()):
      O.vmap[identifier] = "cmn." + prepend_identifier_if_necessary(identifier)
    elif (fdecl.is_save()):
      O.vmap[identifier] = "sve." + prepend_identifier_if_necessary(identifier)
    elif (fdecl.is_intrinsic()):
      if (identifier in ["float", "int", "char"]):
        O.vmap[identifier] = "fem::f" + identifier
      elif (identifier == "iargc"):
        O.vmap[identifier] = "cmn.iargc"
      else:
        O.vmap[identifier] = "fem::" + identifier
    elif (not O.set_vmap_for_callable(identifier=fdecl.id_tok.value)):
      O.set_vmap_force_local(fdecl=fdecl)
      return False
    return True

  def vmapped(O, fdecl):
    identifier = fdecl.id_tok.value
    result = O.vmap.get(identifier)
    if (result is None):
      O.set_vmap_from_fdecl(fdecl=fdecl)
      result = O.vmap[identifier]
    return result

  def vmapped_callable(O, identifier):
    result = O.vmap.get(identifier)
    if (result is None):
      if (not O.set_vmap_for_callable(identifier=identifier)):
        O.vmap[identifier] = prepend_identifier_if_necessary(identifier)
      result = O.vmap[identifier]
    return result

def called_fproc_needs_cmn(conv_info, called_name):
  called_names = conv_info.fproc.externals_passed_by_arg_identifier.get(
    called_name)
  if (called_names is None):
    called_names = [called_name]
  for called_name in called_names:
    sub_fproc = conv_info.fprocs_by_name.get(called_name)
    if (    sub_fproc is not None
        and sub_fproc.needs_cmn
        and not sub_fproc.conv_hook.ignore_common_and_save):
      return True
  return False

def cmn_needs_to_be_inserted(conv_info, prev_tok):
  if (    prev_tok is not None
      and prev_tok.is_identifier()
      and conv_info.fprocs_by_name is not None
      and conv_info.fproc is not None):
    fdecl = conv_info.fproc.get_fdecl(id_tok=prev_tok)
    if (fdecl.is_user_defined_callable()):
      if (called_fproc_needs_cmn(
            conv_info=conv_info,
            called_name=prev_tok.value)):
        return True
  return False

def convert_power(conv_info, tokens):
  fun = "fem::pow"
  pow_tok = tokens[1]
  if (pow_tok.is_integer()):
    if (pow_tok.value == "1"):
      fun = ""
      tokens = tokens[:1]
    elif (pow_tok.value in ["2","3","4"]):
      fun = "fem::pow%s" % pow_tok.value
      tokens = tokens[:1]
  return fun + "(" + convert_tokens(
    conv_info=conv_info, tokens=tokens, commas=True) + ")"

def convert_tokens(conv_info, tokens, commas=False, had_str_concat=None):
  result = []
  rapp = result.append
  prev_tok = None
  prev_prev_tok = None
  if (had_str_concat is None):
    had_str_concat = mutable(value=False)
  from fable.tokenization import group_power
  for tok in group_power(tokens=tokens):
    if (tok.value == "dsign"):
      final = "_MPLAPACK_REPLACE_"
      rapp(final)
    if (tok.value == "dcmplx"):
      final = "_MPLAPACK_REPLACE_"
      rapp(final)
    if (tok.value == "dconjg"):
      final = "_MPLAPACK_REPLACE_"
      rapp(final)
    if (tok.value == "dsqrt"):
      final = "_MPLAPACK_REPLACE_"
      rapp(final)
    if (tok.value == "dble"):
      final = "_MPLAPACK_REPLACE_"
      rapp(final)
    if (tok.value == "dimag"):
      final = "_MPLAPACK_REPLACE_"
      rapp(final)
    if (tok.is_seq()):
      if (    len(tok.value) == 2
          and tok.value[0].is_op_with(value="*")
          and tok.value[1].is_integer()):
        rapp("star /* %s UNHANDLED */" % tok.value[1].value)
      else:
        rapp(convert_tokens(
          conv_info=conv_info,
          tokens=tok.value,
          commas=False,
          had_str_concat=had_str_concat))
    elif (tok.is_parentheses()):
      ed = ")"
      if (cmn_needs_to_be_inserted(conv_info=conv_info, prev_tok=prev_tok)):
        if (len(tok.value) != 0) and (len(tok.value[0].value) != 0):
          op = "(cmn, "
        else:
          op = "(cmn"
      else:
        if (prev_tok is not None and prev_tok.is_identifier() and prev_tok.value != "min" and prev_tok.value != "dsign" and prev_tok.value != "dsqrt" and prev_tok.value != "dble" and prev_tok.value != "dimag" and prev_tok.value != "max" and prev_tok.value != "abs"  and prev_tok.value != "sqrt" and prev_tok.value != "log2" and prev_tok.value != "pow2" and prev_tok.value != "dlamch" and prev_tok.value != "int" and prev_tok.value != "lsame" and prev_tok.value != "lsamen" and prev_tok.value != "ddot" and prev_tok.value != "pow2" and prev_tok.value != "disnan" and prev_tok.value != "mod" and prev_tok.value != "dabs" and prev_tok.value != "cdabs" and prev_tok.value != "dconjg" and prev_tok.value != "dcabs1" and prev_tok.value != "dcmplx" and prev_tok.value != "ceiling" and prev_tok.value != "cos" and prev_tok.value != "sin" and prev_tok.value != "ilaenv"  and prev_tok.value != "disnan" and prev_tok.value != "dlaisnan" and prev_tok.value != "idamax" and prev_tok.value != "ladiv2" and prev_tok.value != "zladiv2" and prev_tok.value != "sign" and prev_tok.value != "iladlr" and prev_tok.value != "iladlc" and prev_tok.value != "stack" and prev_tok.value != "log" and prev_tok.value != "dnrm2" and prev_tok.value != "dlapy2" and prev_tok.value != "dlanst" and prev_tok.value != "dlansy" and prev_tok.value != "idamax" and prev_tok.value != "izamax" and prev_tok.value != "conjg" and prev_tok.value != "dznrm2" and prev_tok.value != "dlapy3"  and prev_tok.value != "atan2" and prev_tok.value != "zdotc" and prev_tok.value != "ilazlr" and prev_tok.value != "zladiv" and prev_tok.value != "zlanhp" and prev_tok.value != "dlansb" and prev_tok.value != "zlanhb" and prev_tok.value != "dlansp" and prev_tok.value != "dlange" and prev_tok.value != "dlantb" and prev_tok.value != "dlantp" and prev_tok.value != "dlantr" and prev_tok.value != "zlanhb" and prev_tok.value != "zlanhe" and prev_tok.value != "zlantb" and prev_tok.value != "zlantp" and prev_tok.value != "zlanhp" and prev_tok.value != "zlanht" and prev_tok.value != "zlansp" and prev_tok.value != "zlansy" and prev_tok.value != "zlange" and prev_tok.value !="dasum" and prev_tok.value !="dla_porcond" and prev_tok.value != "zdotu" and prev_tok.value != "zla_gbrcond_c" and prev_tok.value != "zla_gbrcond_x" and prev_tok.value != "zla_gbrpvgrw" and prev_tok.value != "zla_gercond_c" and prev_tok.value != "zla_gercond_x" and prev_tok.value != "zla_gerpvgrw" and prev_tok.value != "zla_hercond_c" and prev_tok.value != "zla_hercond_x" and prev_tok.value != "zla_herpvgrw" and prev_tok.value != "zla_porcond_c" and prev_tok.value != "zla_porcond_x" and prev_tok.value != "zla_porpvgrw" and prev_tok.value != "zla_syrcond_c" and prev_tok.value != "zla_syrcond_x" and prev_tok.value != "zla_syrpvgrw" and prev_tok.value != "zlangb" and prev_tok.value != "zlangt" and prev_tok.value != "zlanhs" and prev_tok.value != "zlantr" and prev_tok.value != "dzasum" and prev_tok.value != "dzsum1" and prev_tok.value != "dla_gbrcond" and prev_tok.value != "dla_gbrpvgrw" and prev_tok.value != "dla_gercond" and prev_tok.value != "dla_gerpvgrw" and prev_tok.value != "dla_porpvgrw" and prev_tok.value != "dla_syrcond" and prev_tok.value != "dla_syrpvgrw" and prev_tok.value != "dladiv2" and prev_tok.value != "dlaneg" and prev_tok.value != "dlangb" and prev_tok.value != "dlangt" and prev_tok.value != "dlanhs" and prev_tok.value != "rpvgrw" and prev_tok.value != "fchar" and prev_tok.value != "ichar" and prev_tok.value != "cabs1" and prev_tok.value != "cabs2" and prev_tok.value != "abs2" and prev_tok.value != "abs1" and prev_tok.value != "iparam2stage" and prev_tok.value != "ceiling" and prev_tok.value != "exp" and prev_tok.value != "zlarnd" and prev_tok.value != "dlarnd" and prev_tok.value != "dlaran" and prev_tok.value != "nint" and prev_tok.value != "dlog" and prev_tok.value != "ilatrans" and prev_tok.value != "ilauplo" and prev_tok.value != "maxloc" and prev_tok.value != "dlamc3" and prev_tok.value != "log10" and prev_tok.value != "_Rslect" and prev_tok.value != "Rslect" and prev_tok.value != "zslect" and prev_tok.value != "_Cslect" and prev_tok.value != "dslect" and prev_tok.value != "zslect" and prev_tok.value != "zlctes" and prev_tok.value != "dlctes" and prev_tok.value != "clctes"): 
          a=convert_tokens(conv_info=conv_info,tokens=tok.value,commas=True)
          aa = a.split(',')
          if (len(aa)==1):  # one dimensional array
            if(len(aa[0].split(" "))==1):
              op = "["
              ed ="-1]"
            else:
              op = "[("
              ed =")-1]"
            final = op + a.strip() + ed
          else:
            aa0=aa[0].split()
            aa1=aa[1].split()
            if (len(aa)==2): # two dimensional array with i,j
              op = "["
              ed ="]"
              if len(aa0)!=1:
                aaa="((" + aa[0].strip() + ")-1)"
              else:
                aaa="(" + aa[0].strip() + "-1)"
              if len(aa1)!=1:
                bbb="((" + aa[1].strip() + ")-1)*ld" + prev_tok.value
              else:
                bbb="(" + aa[1].strip() + "-1)*ld" + prev_tok.value
              remove0= "+(1-1)*ld" + prev_tok.value
              remove1= "(1-1)+"
              final = op + aaa + "+" + bbb + ed
              _final = final.replace(remove0,'')
              __final = _final.replace(remove1,'')
              rapp(__final)
              continue
            elif (len(aa)==3): # XXX two dimensional array with A[max(i,j),k] , A[max(i, max(j,k)] type
              aa2=aa[2].split()
              op = "["
              ed ="]"
              if '(' in aa[0] and ')' in aa[1]:
                aaaa = aa[0].strip() + ',' + aa[1].strip()
                if len(aa2)!=1:
                  bbbb = '(' + aa[2].strip() + ')'
                else:
                  bbbb = aa[2].strip()
              elif '(' in aa[1] and ')' in aa[2]:
                if len(aa0)!=1:
                  aaaa = '(' + aa[0].strip() + ')'
                else:
                  aaaa = aa[0].strip()
                bbbb = aa[1].strip() + ',' + aa[2].strip()
              else:
                final = aa[0].strip() + ',' + aa[1].strip() + ',' + aa[2].strip()
                rapp(final)
                continue
              final =op + "(" + aaaa + "-1)" + "+" + "(" + bbbb + "-1)*ld" + prev_tok.value + ed
              remove0= "+(1-1)*ld" + prev_tok.value
              remove1= "(1-1)+"
              _final = final.replace(remove0,'')
              __final = _final.replace(remove1,'')
              rapp(__final)
              continue
            else:
              aaa=','.join(aa)
              op = "("
              ed = ")"
              final = op + aaa + ed
              rapp(final)
              continue
        elif (prev_tok is not None and prev_tok.is_identifier() and prev_tok.value == "dble"):
          a=convert_tokens(conv_info=conv_info,tokens=tok.value,commas=True)
          op = "("
          ed = ")"
          ed2 = ".real()"
          aa = delete_brackets(a)
          if aa.isalnum():
            final = a + ed2
          else:
            final = op + a + ed + ed2
          rapp(final)
          continue
        elif (prev_tok is not None and prev_tok.is_identifier() and prev_tok.value == "dimag"):
          a=convert_tokens(conv_info=conv_info,tokens=tok.value,commas=True)
          op = "("
          ed = ")"
          ed2 = ".imag()"
          aa = delete_brackets(a)
          if aa.isalnum():
            final = a + ed2
          else:
            final = op + a + ed + ed2
          rapp(final)
          continue
        elif (prev_tok is not None and prev_tok.is_identifier() and prev_tok.value == "max"):
          a=convert_tokens(conv_info=conv_info,tokens=tok.value,commas=True)
          aa = a.split(',')
          if (len(aa)==2):
            op = "("
            ed = ")"
          else:
            op = "({"
            ed = "})"
          final = op + a + ed
          rapp(final)
          continue
        elif (prev_tok is not None and prev_tok.is_identifier() and prev_tok.value == "min"):
          a=convert_tokens(conv_info=conv_info,tokens=tok.value,commas=True)
          aa = a.split(',')
          if (len(aa)==2):
            op = "("
            ed = ")"
          else:
            op = "({"
            ed = "})"
          final = op + a + ed
          rapp(final)
          continue
        else: # not an array
          op = "("
          ed =")"
      rapp(op + convert_tokens(
        conv_info=conv_info,
        tokens=tok.value,
        commas=True,
        had_str_concat=had_str_concat) + ed )
    elif (tok.is_implied_do()):
      raise AssertionError
    elif (tok.is_power()):
      rapp(convert_power(conv_info=conv_info, tokens=tok.value))
    else:
      rapp(convert_token(
        vmap=conv_info.vmap,
        leading=(len(result)==0),
        tok=tok,
        had_str_concat=had_str_concat))        
    prev_tok = tok
  if (commas):
    return ", ".join(result)
  return "".join(result)

def convert_to_int_literal(tokens):
  assert tokens is not None and len(tokens) != 0
  if (len(tokens) != 1 or not tokens[0].is_integer()):
    tokens[0].raise_not_supported()
  return int(strip_leading_zeros(string=tokens[0].value))

def convert_data_type(conv_info, fdecl, crhs):
  _csize="dummy"
  if (fdecl.data_type is None):
    assert conv_info.fproc is not None
    fdecl.data_type = conv_info.fproc.implicit.get(fdecl.id_tok.value[0])
    if (fdecl.data_type is None):
      raise fdecl.id_tok.raise_semantic_error(msg="Missing data type")
  if (isinstance(fdecl.data_type, str)):
    data_type_code = fdecl.data_type
  else:
    data_type_code = fdecl.data_type.value
  size_tokens = fdecl.size_tokens
  dim_tokens = fdecl.dim_tokens
  if (data_type_code == "character"):
    if (size_tokens is None):
      csize = "1"
    else:
      csize = convert_tokens(conv_info=conv_info, tokens=size_tokens)
    ctype = "char"
    _csize = csize
    if (crhs is None):
      crhs = "fem::char0"
  else:
    def convert_to_ctype_with_size(ctype):
      if (size_tokens is None):
        return ctype
      return "fem::%s_star_%d" % (data_type_code, convert_to_int_literal(
        tokens=size_tokens))
    if (data_type_code == "logical"):
      ctype = convert_to_ctype_with_size(ctype="bool")
    elif (data_type_code == "integer"):
      ctype = convert_to_ctype_with_size(ctype="int")
    elif (data_type_code == "real"):
      ctype = convert_to_ctype_with_size(ctype="float")
    elif (data_type_code == "doubleprecision"):
      if (size_tokens is not None):
        size_tokens[0].raise_syntax_error()
      ctype = "double"
    elif (data_type_code == "byte"):
      if (size_tokens is not None):
        size_tokens[0].raise_syntax_error()
      ctype = "char"
    elif (data_type_code == "complex"):
      if (size_tokens is None):
        ctype = "std::complex<float>"
        if (crhs is None):
          crhs = "fem::float0"
      else:
        sz = convert_to_int_literal(tokens=size_tokens)
        if (sz == 8):
          ctype = "std::complex<float>"
          if (crhs is None):
            crhs = "fem::float0"
          else:
            if crhs[0] == '(' and crhs[-1] == ')':
              crhs = "COMPLEX" + crhs
        elif (sz == 16):
          ctype = "std::complex<double>"
          if (crhs is None):
            crhs = "0.0" #fem::double0
          else:
            if crhs[0] == '(' and crhs[-1] == ')':
              crhs = "COMPLEX" + crhs
        elif (sz == 32):
          ctype = "std::complex<long double>"
          if (crhs is None):
            crhs = "fem::long_double0"
        else:
          size_tokens[0].raise_not_supported()
    elif (data_type_code == "doublecomplex"):
      if (size_tokens is not None):
        size_tokens[0].raise_not_supported()
      ctype = "std::complex<double>"
      if (crhs is None):
        crhs = "0.0" #fem::double0
    else:
      raise RuntimeError(
        "Not implemented: data_type_code = %s" % data_type_code)
  return ctype, crhs, _csize

def convert_dims(conv_info, dim_tokens):
  need_origin = False
  for tokens in dim_tokens:
    for tok in tokens.value:
      if (tok.is_op_with(value=":")):
        need_origin = True
        break
  dims = []
  for i_dim,tokens in enumerate(dim_tokens):
    cdim = convert_tokens(conv_info=conv_info, tokens=tokens.value)
    if (cdim == "star "):
      cdim = "star"
    else:
      cdim = cdim.replace(",  * ", ", star") # XXX
    if (need_origin):
      dims.append("dim%d(%s)" % (i_dim+1, cdim))
    else:
      dims.append(cdim)
  if (need_origin):
    result = ".".join(dims)
  else:
    result = "[" + " * ".join(dims) + "]"
  return result

def parenthesize_if_necessary(expr):
  from fable import unsigned_integer_scan, identifier_scan
  if (   unsigned_integer_scan(code=expr) == len(expr)
      or identifier_scan(code=expr) == len(expr)):
    return expr
  return "(" + expr + ")"

def convert_dim_to_static_size(conv_info, tokens):
  def conv(toks):
    return convert_tokens(conv_info=conv_info, tokens=toks)
  for i,tok in enumerate(tokens):
    if (tok.is_op_with(value=":")):
      f = parenthesize_if_necessary(expr=conv(toks=tokens[:i]))
      l = parenthesize_if_necessary(expr=conv(toks=tokens[i+1:]))
      return "%s-%s+1" % (l, f)
  return "%s" % conv(toks=tokens)

def convert_dims_to_static_size(conv_info, dim_tokens):
  terms = []
  for tokens in dim_tokens:
    terms.append(convert_dim_to_static_size(conv_info, tokens=tokens.value))
  if (len(terms) == 1):
    return terms[0]
  return " * ".join([parenthesize_if_necessary(expr=expr) for expr in terms])

def convert_data_type_and_dims(conv_info, fdecl, crhs, force_arr=False):
  ctype, crhs, _csize = convert_data_type(conv_info=conv_info, fdecl=fdecl, crhs=crhs)
  if _csize == "":
    _csize ="dummy"
  _ctype =""
  dt = fdecl.dim_tokens
  cdims = None
  cfill0 = "fem::fill0"
  if (dt is not None):
    atype = None
    if (    not force_arr
        and conv_info.arr_nd_size_max is not None
        and len(dt) <= 3):
      vals = conv_info.fproc.eval_dimensions_simple(
        dim_tokens=dt, allow_power=False)
      if (vals.count(None) == 0):
        sz = product(vals)
        if (sz <= abs(conv_info.arr_nd_size_max)):
          from fable.read import dimensions_are_simple
          if (dimensions_are_simple(dim_tokens=dt)):
            cdims = Auto
            t = convert_tokens(conv_info=conv_info, tokens=dt, commas=True)
            __t =t.replace(",", "*");
            _csize = __t
          else:
            t = ", ".join(["%d" % v for v in vals])
            __t = "* ".join(["%d" % v for v in vals])
          _t = "%s" % __t
          _ctype = "%s" % ctype
          _csize = t.replace(",", "*");
          t = "%s, %s" % (t, ctype)
          if (t.endswith(">")): templs = " "
          else:                 templs = ""
          atype = "arr_%dd<%s%s>" % (len(dt), t, templs)
          if (conv_info.arr_nd_size_max < 0):
            cfill0 = "fem::no_fill0"
    if (atype is None):
      if (len(dt) != 1):
        t = "%s, %d" % (ctype, len(dt))
      else:
        t = ctype
      if (t.endswith(">")): templs = " "
      else:                 templs = ""
      atype = "arr<%s%s>" % (t, templs)
      _ctype = ctype
    ctype = atype
    if (cdims is None):
      cdims = convert_dims(conv_info=conv_info, dim_tokens=dt)
  return ctype, cdims, crhs, cfill0, _csize, _ctype

def ad_hoc_change_arr_to_arr_ref(ctype, cconst=""):
  return ctype.replace("arr<", "arr_%sref<" % cconst, 1)

def zero_shortcut_if_possible(ctype):
  if (ctype.startswith("fem::")):
    if (ctype.endswith(">")): s = " "
    else:                     s = ""
    return "fem::zero<%s%s>()" % (ctype, s)
#  return "%s" % ctype
  return "0"

def convert_declaration(rapp, conv_info, fdecl, crhs, const):
  ctype, cdims, crhs, cfill0, _csize, _ctype  = convert_data_type_and_dims(
    conv_info=conv_info, fdecl=fdecl, crhs=crhs)
  vname = conv_info.vmapped(fdecl=fdecl)
  if (cdims is None):
    if (crhs is None): crhs = zero_shortcut_if_possible(ctype=ctype)
    def const_qualifier():
      if (const): return "const "
      return ""
    if _csize == "" or _csize == "dummy":
      if ctype == "double":
          rapp("%s%s %s = %s;" % (const_qualifier(), "REAL", vname, crhs))
      elif ctype == "int":
          rapp("%s%s %s = %s;" % (const_qualifier(), "INTEGER", vname, crhs))
      else:
          rapp("%s%s %s = %s;" % (const_qualifier(), ctype, vname, crhs))
    else:
      if _csize =="1":
        rapp("%s%s %s;" % (const_qualifier(), ctype, vname))
      else:
        rapp("%s%s %s[%s] = %s;" % (const_qualifier(), ctype, vname, _csize, crhs))
    return False
  if (cdims is Auto):
    if _ctype == "":
      rapp("%s %s;" % (ctype, vname))
    else:
      if _ctype == "char" and _csize =="1":
        rapp("%s %s;" % (_ctype, vname ))
      else:
        rapp("%s %s [%s];" % (_ctype, vname, _csize ))
  else:
    if _ctype == "":
      rapp("%s %s %s;" % (ctype, vname, cdims))
    else:
      rapp("%s %s %s;" % (_ctype, vname, cdims))
  return True

class scope(object):

  __slots__ = [
    "parent",
    "opening_text",
    "auto_close_parent",
    "closing_text",
    "data",
    "insert_point",
    "trailing_statement_label_index",
    "tail"]

  def __init__(O, parent, opening_text=None, auto_close_parent=False):
    O.parent = parent
    O.opening_text = opening_text
    O.auto_close_parent = auto_close_parent
    O.closing_text = None
    O.data = []
    O.insert_point = None
    O.trailing_statement_label_index = None
    O.tail = None

  def current_point(O):
    return len(O.data)

  def point_is_current(O, point):
    return (point == len(O.data))

  def remember_insert_point(O):
    O.insert_point = len(O.data)

  def insert_point_is_current(O):
    return O.point_is_current(point=O.insert_point)

  def top_append(O, obj):
    assert O.insert_point is not None
    O.data.insert(O.insert_point, obj)
    O.insert_point += 1

  def append(O, obj):
    O.data.append(obj)
    O.trailing_statement_label_index = None

  def append_statement_label(O, label):
    O.trailing_statement_label_index = len(O.data)
    O.data.append("statement_%s:" % label)

  def append_comment(O, line):
    O.data.append(line)

  def open_nested_scope(O, opening_text, auto_close_parent=False):
    O.trailing_statement_label_index = None
    return scope(
      parent=O, opening_text=opening_text, auto_close_parent=auto_close_parent)

  def finalize(O):
    if (O.trailing_statement_label_index is not None):
      O.data[O.trailing_statement_label_index] += ";"

  def close_nested_scope(O):
    assert O.opening_text is not None
    assert O.closing_text is None
    assert O.tail is None
    O.finalize()
    O.closing_text = ["}"]
    head = O
    while (head.parent.tail is head):
      head = head.parent
    head.parent.data.append(head)
    if (head.auto_close_parent):
      return head.parent.close_nested_scope()
    return head.parent

  def attach_tail(O, opening_text):
    assert O.opening_text is not None
    assert O.closing_text is None
    assert O.tail is None
    O.finalize()
    O.closing_text = ["}"]
    O.tail = scope(parent=O, opening_text=opening_text)
    return O.tail

  def collect_text(O, callback, indent="  "):
    for obj in O.data:
      if (isinstance(obj, scope)):
        curr = obj
        while (curr is not None):
          for text in curr.opening_text:
            callback(indent+text)
          curr.collect_text(callback=callback, indent=indent+"  ")
          for text in curr.closing_text:
            callback(indent+text)
          curr = curr.tail
      else:
        callback(indent+obj)

def convert_io_statement_with_err(
      conv_info,
      curr_scope,
      io_function,
      io_function_specialization,
      io_call_args,
      iolist):
  if (iolist.err is None):
    io_scope = curr_scope
  else:
    io_scope = curr_scope.open_nested_scope(opening_text=["try {"])
  io_scope.append("cmn.io.%s%s(%s)" % (
    io_function, io_function_specialization, io_call_args))
  for slot in iolist.chain:
    tokens = getattr(iolist, slot)
    if (tokens is not None):
      carg = convert_tokens(conv_info=conv_info, tokens=tokens)
      io_scope.append("  .%s(%s)" % (slot, carg))
  io_scope.data[-1] += ";"
  if (io_scope is not curr_scope):
    io_scope.close_nested_scope()
    catch_scope = curr_scope.open_nested_scope(
      opening_text=["catch (fem::io_err const&) {"])
    clabel = convert_tokens(conv_info=conv_info, tokens=iolist.err)
    catch_scope.append("goto statement_%s;" % clabel)
    catch_scope.close_nested_scope()

def is_simple_do_last(tokens):
  i = 0
  if (len(tokens) == 2 and tokens[0].is_unary_plus_or_minus()):
    i = 1
  if (i+1 == len(tokens)):
    tok = tokens[i]
    return tok.is_identifier() or tok.is_integer()
  return False

def convert_to_fem_do(conv_info, parent_scope, i_tok, fls_tokens):
  assert 2 <= len(fls_tokens) <= 3
  i = convert_token(vmap=conv_info.vmap, leading=True, tok=i_tok)
  f = convert_tokens(conv_info=conv_info, tokens=fls_tokens[0].value)
  l = convert_tokens(conv_info=conv_info, tokens=fls_tokens[1].value)
  if (len(fls_tokens) == 3):
    s = convert_tokens(conv_info=conv_info, tokens=fls_tokens[2].value)
    if (s.lstrip('+-').isdigit()):
      if int(s) >=0:
        return parent_scope.open_nested_scope(
         opening_text=["for(%s=%s; %s<=%s; %s=%s+%s) {" % (i, f, i, l, i, i, s)])
      else:
        return parent_scope.open_nested_scope(
         opening_text=["for(%s=%s; %s>=%s; %s=%s%s) {" % (i, f, i, l, i, i, s)])
    else:
      if '-' in s:
        return parent_scope.open_nested_scope(
               opening_text=["for(%s=%s; %s>=%s; %s=%s%s) {" % (i, f, i, l, i, i, s)])
      else:
        return parent_scope.open_nested_scope(
               opening_text=["for(%s=%s; %s<=%s; %s=%s+%s) {" % (i, f, i, l, i, i, s)])
  if (conv_info.fem_do_safe):
    return parent_scope.open_nested_scope(
      opening_text=["for(%s=%s; %s<=%s; %s=%s+1) {" % (i, f, i, l, i, i)])
  if (is_simple_do_last(tokens=fls_tokens[1].value)):
    return parent_scope.open_nested_scope(
      opening_text=["FEM_DO(%s, %s, %s) {" % (i, f, l)])
  scope_for_last = parent_scope.open_nested_scope(opening_text=["{"])
  scope_for_last.append("int fem_do_last = %s;" % l)
  return scope_for_last.open_nested_scope(
    opening_text=["FEM_DO(%s, %s, fem_do_last) {" % (i, f)],
    auto_close_parent=True)

def find_implied_dos(result, tokens):
  assert isinstance(tokens, list)
  for i,tok in enumerate(tokens):
    if (tok.is_seq_or_parentheses()):
      find_implied_dos(result=result, tokens=tok.value)
    elif (tok.is_implied_do()):
      result.append(tok)

def convert_io_loop(
      io_scope, io_op, conv_info, tokens, cbuf=None, had_str_concat=None):
  class cbuffer(object):
    __slots__ = ["strings", "leading"]
    def __init__(O):
      O.strings = []
      O.leading = True
    def append(O, string):
      O.strings.append(string)
      O.leading = False
    def append_comma(O):
      if (len(O.strings) != 0 and O.strings[-1] != ", "):
        O.strings.append(", ")
      O.leading = True
    def remove_trailing_comma(O):
      if (len(O.strings) != 0):
        assert O.strings[-1] == ", "
        O.strings.pop()
        assert O.leading
        O.leading = False
    def append_opening_parenthesis(O):
      O.strings.append("(")
      O.leading = True
    def append_closing_parenthesis(O):
      assert len(O.strings) != 0
      if (O.strings[-1] == ", "):
        O.remove_trailing_comma()
      O.strings.append(")")
    def flush(O):
      O.remove_trailing_comma()
      if (len(O.strings) != 0):
        io_scope.append("%s, %s;" % (io_op, "".join(O.strings)))
        O.strings = []
  if (cbuf is None):
    cbuf = cbuffer()
    owning_cbuf = True
  else:
    owning_cbuf = False
  prev_tok = None
  if (had_str_concat is None):
    had_str_concat = mutable(value=False)
  from fable.tokenization import group_power
  for tok in group_power(tokens=tokens):
    if (tok.is_seq()):
      convert_io_loop(
        io_scope,
        io_op,
        conv_info,
        tokens=tok.value,
        cbuf=cbuf,
        had_str_concat=had_str_concat)
      cbuf.append_comma()
    elif (tok.is_parentheses()):
      cbuf.append_opening_parenthesis()
      if (cmn_needs_to_be_inserted(conv_info=conv_info, prev_tok=prev_tok)):
        cbuf.append("cmn")
        if (len(tok.value) != 0):
          cbuf.append_comma()
      convert_io_loop(
        io_scope,
        io_op,
        conv_info,
        tokens=tok.value,
        cbuf=cbuf,
        had_str_concat=had_str_concat)
      cbuf.append_closing_parenthesis()
    elif (tok.is_implied_do()):
      cbuf.flush()
      from fable.tokenization import implied_do_info
      idi = implied_do_info(tokens=tok.value)
      do_scope = convert_to_fem_do(
        conv_info=conv_info,
        parent_scope=io_scope,
        i_tok=idi.id_tok,
        fls_tokens=idi.fls_tokens)
      convert_io_loop(
        io_scope=do_scope,
        io_op=io_op,
        conv_info=conv_info,
        tokens=tok.value[:idi.dlist_size],
        had_str_concat=had_str_concat)
      do_scope.close_nested_scope()
      return
    elif (tok.is_power()):
      cbuf.append(convert_power(conv_info=conv_info, tokens=tok.value))
    else:
      cbuf.append(convert_token(
        vmap=conv_info.vmap,
        leading=cbuf.leading,
        tok=tok,
        had_str_concat=had_str_concat))
    prev_tok = tok
  if (owning_cbuf):
    cbuf.flush()

def equivalence_align_with_arg(conv_info, top_scope, identifier, tok_seq):
  assert tok_seq.is_seq()
  tokens = tok_seq.value
  assert len(tokens) > 0
  if (len(tokens) == 1):
    return ""
  cindices = []
  for i in range(1,len(tokens)):
    tok = tokens[i]
    if (i == 3 or not tok.is_parentheses()):
      tok.raise_semantic_error()
    declare_identifiers_parameter_recursion(
      conv_info=conv_info,
      top_scope=top_scope,
      curr_scope=top_scope,
      tokens=tokens[i].value)
    cindices.append(convert_tokens(
      conv_info=conv_info, tokens=tokens[i].value, commas=True))
  fdecl = conv_info.fproc.fdecl_by_identifier[identifier]
  if (len(cindices) == 1):
    if (fdecl.dim_tokens is not None):
      return "arr_index(%s)" % cindices[0]
    if (fdecl.data_type.value != "character"):
      tok_seq.raise_semantic_error()
    return "str_index(%s)" % cindices[0]
  if (   len(cindices) != 2
      or fdecl.data_type.value != "character"
      or fdecl.dim_tokens is None):
    tok_seq.raise_semantic_error()
  return "arr_index(%s)(%s)" % tuple(cindices)

def cconst(fdecl, short):
  if (fdecl.is_modified): return ""
  if (short): return "c"
  return " const"

def convert_to_mbr_bind(
      conv_info,
      top_scope,
      variant_bind_chain,
      mbr_buffer,
      bind_buffer,
      identifier):
  fdecl = conv_info.fproc.fdecl_by_identifier[identifier]
  ctype = convert_data_type(conv_info=conv_info, fdecl=fdecl, crhs=None)[0]
  if (fdecl.dim_tokens is None):
    cdims_parens = ""
  else:
    declare_identifiers_parameter_recursion(
      conv_info=conv_info,
      top_scope=top_scope,
      curr_scope=top_scope,
      tokens=fdecl.dim_tokens)
    cdims = convert_dims(conv_info=conv_info, dim_tokens=fdecl.dim_tokens)
    cdims_parens = "(" + cdims + ")"
  identifier = fdecl.id_tok.value
  ctype_targ = ctype
  if (ctype_targ.endswith(">")):
    ctype_targ += " "
  mbr_buffer.append("mbr<%s> %s%s;" % (ctype_targ, identifier, cdims_parens))
  conv_info.set_vmap_force_local(fdecl=fdecl)
  vname = conv_info.vmapped(fdecl=fdecl)
  if (fdecl.use_count == 0):
    pr = "/* "
    eq = "*/"
    clm = " */ "
    prm = " /* "
  else:
    pr = ""
    eq = "="
    clm = ""
    prm = ""
  if (fdecl.dim_tokens is None):
    if (fdecl.data_type.value == "character"):
      binding = "%sstr_%sref %s %s %s.bind_str();" % (
        pr, cconst(fdecl=fdecl, short=True), vname, eq, variant_bind_chain)
    else:
      binding = "%s%s%s& %s %s %s.bind<%s>();" % (
        pr, ctype, cconst(fdecl=fdecl, short=False), vname,
        eq, variant_bind_chain, ctype_targ)
  else:
    if (fdecl.data_type.value == "character"):
      ref_dim = "%d" % len(fdecl.dim_tokens)
      if (ref_dim == "1"): ref_dim = ""
      binding = "%sstr_arr_%sref<%s> %s(%s%s.bind_str()%s, %s)%s;" % (
        pr, cconst(fdecl=fdecl, short=True), ref_dim, vname,
        clm, variant_bind_chain, prm, cdims, clm)
    else:
      ref_dim = ", %d" % len(fdecl.dim_tokens)
      if (ref_dim == ", 1"): ref_dim = ""
      binding = "%sarr_%sref<%s%s> %s(%s%s.bind<%s>()%s, %s)%s;" % (
        pr, cconst(fdecl=fdecl, short=True), ctype, ref_dim, vname,
        clm, variant_bind_chain, ctype, prm, cdims, clm)
  bind_buffer.append(binding)

def assemble_allocate_line_lists(
      conv_info,
      top_scope,
      variant_bind_chain,
      mbr_buffer,
      bind_buffer,
      allocate_line_lists,
      equiv_tok_cluster,
      identifier):
  if (allocate_line_lists[-1] == [" "]):
    allocate_line_lists.pop()
  i_mbr_by_identifer = {identifier: 0}
  eq_identifiers = [identifier]
  i_block = len(allocate_line_lists)
  allocate_line_lists.append(None)
  for equiv_tok in equiv_tok_cluster:
    align_with = ".align"
    for tok_seq in equiv_tok.value:
      eq_identifier = tok_seq.value[0].value
      i_mbr = i_mbr_by_identifer.get(eq_identifier)
      if (i_mbr is None):
        i_mbr_by_identifer[eq_identifier] = i_mbr = len(eq_identifiers)
        eq_identifiers.append(eq_identifier)
        if (bind_buffer is not None):
          convert_to_mbr_bind(
            conv_info=conv_info,
            top_scope=top_scope,
            variant_bind_chain=variant_bind_chain,
            mbr_buffer=mbr_buffer,
            bind_buffer=bind_buffer,
            identifier=eq_identifier)
      allocate_line_lists.append([
        "    %s<%d>(%s)" % (
          align_with,
          i_mbr+1,
          equivalence_align_with_arg(
            conv_info=conv_info,
            top_scope=top_scope,
            identifier=eq_identifier,
            tok_seq=tok_seq))])
      align_with = " .with"
  allocate_line_lists[i_block] = \
    ["  equivalence(%s)" % (", ".join(eq_identifiers))]
  allocate_line_lists[-1][-1] += ","
  allocate_line_lists.append([" "])

def add_allocate_lines_to_mbr_scope(allocate_line_lists, mbr_buffer):
  if (allocate_line_lists[-1] == [" "]):
    allocate_line_lists.pop()
  if (allocate_line_lists[-1][-1][-1] == ","):
    allocate_line_lists[-1][-1] = allocate_line_lists[-1][-1][:-1]
  if (len(allocate_line_lists) == 1):
    allocate_line_lists[-1][-1] += ";"
  else:
    allocate_line_lists.append([";"])
  for line_list in allocate_line_lists:
    mbr_buffer.append(" ".join(line_list))

def convert_variant_allocate_and_bindings(conv_info, top_scope):
  result_buffers = group_args(
    first_time=[],
    loc_equivalences=[],
    bindings=[])
  equiv_info = conv_info.fproc.equivalence_info()
  equiv_tok_clusters = equiv_info.equiv_tok_clusters
  for common_name,common_fdecl_list in conv_info.fproc.common.items():
    vcn = conv_info.fproc.conv_hook.variant_common_names
    if (vcn is None or common_name not in vcn):
      continue
    top_scope.append(
      "common_variant %s(cmn.common_%s, sve.%s_bindings);" % (
        (common_name,)*3))
    mbr_buffer = []
    result_buffers.first_time.append(mbr_buffer)
    allocate_line_lists = [["%s.allocate()," % common_name]]
    for fdecl in common_fdecl_list:
      identifier = fdecl.id_tok.value
      convert_to_mbr_bind(
        conv_info=conv_info,
        top_scope=top_scope,
        variant_bind_chain=common_name,
        mbr_buffer=mbr_buffer,
        bind_buffer=result_buffers.bindings,
        identifier=identifier)
      equiv_tok_cluster = equiv_info.equiv_tok_cluster_by_identifier.get(
        identifier)
      if (equiv_tok_cluster is None):
        allocate_line_lists[-1].append(identifier+",")
      else:
        assemble_allocate_line_lists(
          conv_info=conv_info,
          top_scope=top_scope,
          variant_bind_chain=common_name,
          mbr_buffer=mbr_buffer,
          bind_buffer=result_buffers.bindings,
          allocate_line_lists=allocate_line_lists,
          equiv_tok_cluster=equiv_tok_cluster,
          identifier=identifier)
    add_allocate_lines_to_mbr_scope(
      allocate_line_lists=allocate_line_lists, mbr_buffer=mbr_buffer)
  #
  cei = conv_info.fproc.classified_equivalence_info()
  if (len(cei.save.equiv_tok_clusters) != 0):
    top_scope.append(
      "save_equivalences sve_equivalences(sve.save_equivalences);")
    mbr_buffer = []
    result_buffers.first_time.append(mbr_buffer)
    mbr_defined_already = set()
    for equiv_tok_cluster in cei.save.equiv_tok_clusters:
      for equiv_tok in equiv_tok_cluster:
        for tok_seq in equiv_tok.value:
          identifier = tok_seq.value[0].value
          if (identifier in mbr_defined_already):
            continue
          mbr_defined_already.add(identifier)
          convert_to_mbr_bind(
            conv_info=conv_info,
            top_scope=top_scope,
            variant_bind_chain="sve_equivalences",
            mbr_buffer=mbr_buffer,
            bind_buffer=result_buffers.bindings,
            identifier=identifier)
    allocate_line_lists = [["sve_equivalences.allocate(),"]]
    for equiv_tok_cluster in cei.save.equiv_tok_clusters:
      assemble_allocate_line_lists(
        conv_info=conv_info,
        top_scope=top_scope,
        variant_bind_chain=None,
        mbr_buffer=mbr_buffer,
        bind_buffer=None,
        allocate_line_lists=allocate_line_lists,
        equiv_tok_cluster=equiv_tok_cluster,
        identifier=equiv_tok_cluster[0].value[0].value[0].value)
    add_allocate_lines_to_mbr_scope(
        allocate_line_lists=allocate_line_lists, mbr_buffer=mbr_buffer)
  #
  if (len(cei.local.equiv_tok_clusters) != 0):
    mbr_defined_already = set()
    for equiv_tok_cluster in cei.local.equiv_tok_clusters:
      for equiv_tok in equiv_tok_cluster:
        for tok_seq in equiv_tok.value:
          identifier = tok_seq.value[0].value
          if (identifier in mbr_defined_already):
            continue
          mbr_defined_already.add(identifier)
          convert_to_mbr_bind(
            conv_info=conv_info,
            top_scope=top_scope,
            variant_bind_chain="loc_equivalences",
            mbr_buffer=result_buffers.loc_equivalences,
            bind_buffer=result_buffers.bindings,
            identifier=identifier)
    allocate_line_lists = [["loc_equivalences.allocate(),"]]
    for equiv_tok_cluster in cei.local.equiv_tok_clusters:
      assemble_allocate_line_lists(
        conv_info=conv_info,
        top_scope=top_scope,
        variant_bind_chain=None,
        mbr_buffer=result_buffers.loc_equivalences,
        bind_buffer=None,
        allocate_line_lists=allocate_line_lists,
        equiv_tok_cluster=equiv_tok_cluster,
        identifier=equiv_tok_cluster[0].value[0].value[0].value)
    add_allocate_lines_to_mbr_scope(
      allocate_line_lists=allocate_line_lists,
      mbr_buffer=result_buffers.loc_equivalences)
  #
  return result_buffers

def convert_data(conv_info, data_init_scope):
  for nlist,clist in conv_info.fproc.data:
    ccs = []
    have_repetitions = False
    tok_types = set()
    for repetition_tok,ctoks in clist:
      i = 0
      if (ctoks[0].is_unary_plus_or_minus() and len(ctoks) > 1):
        i = 1
      tok_types.add(ctoks[i].type())
      cc = convert_tokens(conv_info=conv_info, tokens=ctoks)
      if (repetition_tok is not None):
        have_repetitions = True
        cr = convert_tokens(conv_info=conv_info, tokens=[repetition_tok])
        cc = "%s*datum(%s)" % (cr, cc)
      ccs.append(cc)
    homogeneous_ctype = None
    if (    conv_info.data_specializations
        and not have_repetitions
        and len(tok_types) == 1):
      homogeneous_ctype = {
        "integer": "int",
        "hexadecimal": "int",
        "real": "float",
        "double_precision": "double",
        "logical": "bool",
        "string": "char*",
        "complex": None # TODO
      }.get(list(tok_types)[0])
    def data_values_blocked():
      data_scope.append("fem::data_values data;")
      for i_block in range(0, len(ccs), conv_info.data_values_block_size):
        data_scope.append(
          "data.values, %s;"
            % ", ".join(ccs[i_block:i_block+conv_info.data_values_block_size]))
    def values_for_data_of_type():
      data_scope.append("static const %s values[] = {" % homogeneous_ctype)
      data_scope.append("  %s" % ", ".join(ccs))
      data_scope.append("};")
      if (homogeneous_ctype == "char*"): s = "_str"
      else:                              s = "<%s>" % homogeneous_ctype
      return s
    def have_no_array_targets():
      for tok_seq in nlist:
        if (len(tok_seq.value) == 1):
          fdecl = conv_info.fproc.get_fdecl(id_tok=tok_seq.value[0])
          if (fdecl.dim_tokens is not None):
            return False
      return True
    implied_dos = []
    find_implied_dos(result=implied_dos, tokens=nlist)
    if (len(implied_dos) == 0):
      if (    conv_info.data_specializations
          and len(nlist) == len(ccs)
          and have_no_array_targets()):
        for tok_seq,cc in zip(nlist, ccs):
          cn = convert_tokens(conv_info=conv_info, tokens=tok_seq.value)
          data_init_scope.append("%s = %s;" % (cn, cc))
      else:
        cn = convert_tokens(conv_info=conv_info, tokens=nlist, commas=True)
        if (homogeneous_ctype is None):
          if (len(ccs) <= conv_info.data_values_block_size):
            data_init_scope.append(
              "fem::data((values, %s)), %s;" % (", ".join(ccs), cn))
          else:
            data_scope = data_init_scope.open_nested_scope(opening_text=["{"])
            data_values_blocked()
            data_scope.append("data, %s;" % cn)
            data_scope.close_nested_scope()
        elif (len(ccs) != 1):
          if (    len(conv_info.fproc.data) == 1
              and len(conv_info.fproc.conv_hook.variant_common_names) == 0):
            data_scope = data_init_scope
          else:
            data_scope = data_init_scope.open_nested_scope(opening_text=["{"])
          s = values_for_data_of_type()
          data_scope.append("fem::data_of_type%s(FEM_VALUES_AND_SIZE)," % s)
          data_scope.append("  %s;" % cn)
          if (data_scope is not data_init_scope):
            data_scope.close_nested_scope()
        else:
          data_init_scope.append("%s = %s;" % (cn, ccs[0]))
    else:
      if (    len(conv_info.fproc.data) == 1
          and len(conv_info.fproc.conv_hook.variant_common_names) == 0
          and len(ccs) <= conv_info.data_values_block_size):
        data_scope = data_init_scope
      else:
        data_scope = data_init_scope.open_nested_scope(opening_text=["{"])
      if (homogeneous_ctype is None):
        if (len(ccs) <= conv_info.data_values_block_size):
          data_scope.append(
            "fem::data_values data((values, %s));" % ", ".join(ccs))
        else:
          data_values_blocked()
      else:
        s = values_for_data_of_type()
        data_scope.append("fem::data_of_type%s data(FEM_VALUES_AND_SIZE);" % s)
      convert_io_loop(
        io_scope=data_scope,
        io_op="data",
        conv_info=conv_info,
        tokens=nlist)
      if (data_scope is not data_init_scope):
        data_scope.close_nested_scope()

def declare_identifiers_parameter_recursion(
      conv_info, top_scope, curr_scope, tokens):
  from fable.tokenization import extract_identifiers
  for id_tok in extract_identifiers(tokens=tokens):
    if (id_tok.value in conv_info.vmap):
      continue
    conv_info.vmap[id_tok.value] = None
    fdecl = conv_info.fproc.get_fdecl(id_tok=id_tok)
    declare_identifiers_parameter_recursion(
      conv_info=conv_info, top_scope=top_scope, curr_scope=curr_scope,
      tokens=fdecl.required_parameter_assignment_tokens())
    declare_identifier(
      conv_info=conv_info,
      top_scope=top_scope,
      curr_scope=curr_scope,
      id_tok=id_tok)

def declare_size_dim_identifiers(conv_info, top_scope, curr_scope, fdecl):
  for sd_tokens in [fdecl.size_tokens, fdecl.dim_tokens]:
    if (sd_tokens is None):
      continue
    from fable.tokenization import extract_identifiers
    sd_id_tokens = extract_identifiers(tokens=sd_tokens)
    for sd_id_tok in sd_id_tokens:
      if (sd_id_tok.value == fdecl.id_tok.value):
        sd_id_tok.raise_semantic_error(msg="Recursion in declaration")
      if (sd_id_tok.value in conv_info.vmap):
        continue
      sd_fdecl = conv_info.fproc.get_fdecl(id_tok=sd_id_tok)
      if (sd_fdecl.parameter_assignment_tokens is None):
        sd_crhs = None
      else:
        declare_identifiers_parameter_recursion(
          conv_info=conv_info, top_scope=top_scope, curr_scope=curr_scope,
          tokens=sd_fdecl.parameter_assignment_tokens)
        if (sd_id_tok.value in conv_info.fproc.dynamic_parameters):
          sd_crhs = "cmn.dynamic_params." + sd_id_tok.value
        else:
          sd_crhs = convert_tokens(
            conv_info=conv_info, tokens=sd_fdecl.parameter_assignment_tokens)
      if (not conv_info.set_vmap_from_fdecl(fdecl=sd_fdecl)):
        have_goto = (len(conv_info.fproc.target_statement_labels()) != 0)
        if (have_goto):
          rapp = top_scope.top_append
        else:
          rapp = top_scope.append
        convert_declaration(
          rapp=rapp,
          conv_info=conv_info,
          fdecl=sd_fdecl,
          crhs=sd_crhs,
          const=True)

def simple_equivalence(
      conv_info,
      top_scope,
      curr_scope,
      target_fdecl,
      equiv_tok_cluster):
  assert len(equiv_tok_cluster) != 0
  for equiv_tok in equiv_tok_cluster:
    target_tok_seq = None
    source_tok_seq = None
    source_fdecl = None
    for tok_seq in equiv_tok.value:
      identifier = tok_seq.value[0].value
      if (identifier == target_fdecl.id_tok.value):
        assert target_tok_seq is None
        target_tok_seq = tok_seq
      else:
        fdecl = conv_info.fproc.get_fdecl(id_tok=tok_seq.value[0])
        if (fdecl is not None and fdecl.is_common()):
          assert source_tok_seq is None
          source_tok_seq = tok_seq
          source_fdecl = fdecl
    if (    target_tok_seq is not None
        and source_tok_seq is not None):
      break
  else:
    raise AssertionError
  conv_info.set_vmap_force_local(fdecl=target_fdecl)
  if (conv_info.vmap.get(source_fdecl.id_tok.value) is None):
    declare_identifier(
      conv_info=conv_info,
      top_scope=top_scope,
      curr_scope=curr_scope,
      id_tok=source_fdecl.id_tok)
  crhs = convert_tokens(conv_info=conv_info, tokens=[source_tok_seq])
  se = "// SIMPLE EQUIVALENCE"
  if (target_fdecl.data_type.value == "character"):
    clen = convert_tokens(
      conv_info=conv_info, tokens=target_fdecl.size_tokens)
    if (target_fdecl.dim_tokens is None):
      return "str_%sref %s(%s, %s); %s" % (
        cconst(fdecl=target_fdecl, short=True),
        target_fdecl.id_tok.value, crhs, clen, se)
    cdims = convert_dims(
      conv_info=conv_info, dim_tokens=target_fdecl.dim_tokens)
    return "str_arr_%sref<%d> %s(%s, %s, %s); %s" % (
      cconst(fdecl=target_fdecl, short=True),
      len(target_fdecl.dim_tokens),
      target_fdecl.id_tok.value, crhs, clen, cdims, se)
  ctype, cdims = convert_data_type_and_dims(
    conv_info=conv_info, fdecl=target_fdecl, crhs=None, force_arr=True)[:2]
  if (cdims is None):
    return "%s%s& %s = %s; %s" % (
      ctype, cconst(fdecl=target_fdecl, short=False),
      target_fdecl.id_tok.value, crhs, se)
  return "%s %s(%s, %s); %s" % (
    ad_hoc_change_arr_to_arr_ref(
      ctype=ctype, cconst=cconst(fdecl=target_fdecl, short=True)),
    target_fdecl.id_tok.value, crhs, cdims, se)

def declare_identifier(conv_info, top_scope, curr_scope, id_tok, crhs=None):
  fdecl = conv_info.fproc.get_fdecl(id_tok=id_tok)
  conv_info.set_vmap_from_fdecl(fdecl=fdecl)
  have_goto = (len(conv_info.fproc.target_statement_labels()) != 0)
  def get_rapp():
    if (have_goto):
      return top_scope.top_append
    return top_scope.append
  if (not fdecl.is_common()):
    equiv_tok_cluster = conv_info.fproc.equivalence_info() \
      .equiv_tok_cluster_by_identifier.get(id_tok.value)
    if (equiv_tok_cluster is not None):
      rapp = get_rapp()
      rapp(simple_equivalence(
        conv_info=conv_info,
        top_scope=top_scope,
        curr_scope=curr_scope,
        target_fdecl=fdecl,
        equiv_tok_cluster=equiv_tok_cluster))
      return crhs is not None
  if (fdecl is not None
        and (   fdecl.is_local()
             or fdecl.is_parameter())):
    const = False
    have_crhs = (crhs is not None)
    if (have_goto or curr_scope != top_scope):
      crhs = None
    if (crhs is None):
      if (fdecl.parameter_assignment_tokens is not None):
        declare_identifiers_parameter_recursion(
          conv_info=conv_info, top_scope=top_scope, curr_scope=curr_scope,
          tokens=fdecl.parameter_assignment_tokens)
        if (id_tok.value in conv_info.fproc.dynamic_parameters):
          crhs = "cmn.dynamic_params." + prepend_identifier_if_necessary(
            id_tok.value)
        else:
          crhs = convert_tokens(
            conv_info=conv_info, tokens=fdecl.parameter_assignment_tokens)
        const = True
    elif (fdecl.parameter_assignment_tokens is not None):
      id_tok.raise_semantic_error(
        msg="Assignment to PARAMETER %s" % id_tok.value)
    rapp = get_rapp()
    declare_size_dim_identifiers(
      conv_info=conv_info, top_scope=top_scope, curr_scope=curr_scope,
      fdecl=fdecl)
    result = convert_declaration(
      rapp=rapp,
      conv_info=conv_info,
      fdecl=fdecl,
      crhs=crhs,
      const=const)
    if (have_crhs and (have_goto or curr_scope != top_scope)):
      result = True
    return result
  identifier = id_tok.value
  def get_common_name_if_cast_is_needed():
    if (not fdecl.is_common()): return None
    if (conv_info.converted_commons_info is None): return None
    common_names = conv_info.converted_commons_info.member_registry.get(
      identifier)
    if (common_names is None): return None
    if (len(common_names) < 2): return None
    return conv_info.fproc.common_name_by_identifier().get(identifier)
  common_name = get_common_name_if_cast_is_needed()
  if (common_name is not None):
    src_var = "static_cast<common_%s&>(cmn).%s" % (
      common_name, prepend_identifier_if_necessary(identifier))
  else:
    src_var = conv_info.vmap[identifier]
  if (fdecl.dim_tokens is not None):
    conv_info.vmap[identifier] = prepend_identifier_if_necessary(identifier)
    if (fdecl.data_type.value == "character"):
      ctype = "str_arr_%sref<%d>" % (
        cconst(fdecl=fdecl, short=True), len(fdecl.dim_tokens))
    else:
      ctype = convert_data_type_and_dims(
        conv_info=conv_info, fdecl=fdecl, crhs=None, force_arr=True)[0]
      ctype = ad_hoc_change_arr_to_arr_ref(
        ctype=ctype, cconst=cconst(fdecl=fdecl, short=True))
    declare_size_dim_identifiers(
      conv_info=conv_info,
      top_scope=top_scope,
      curr_scope=top_scope,
      fdecl=fdecl)
    cdims = convert_dims(conv_info=conv_info, dim_tokens=fdecl.dim_tokens)
    if (common_name is not None):
      src_var = "static_cast<common_%s&>(cmn).%s" % (
        common_name, prepend_identifier_if_necessary(identifier))
    rapp = get_rapp()
    rapp("%s %s(%s, %s);" % (
      ctype, prepend_identifier_if_necessary(identifier), src_var, cdims))
  elif (   common_name is not None
        or (fdecl.use_count > 1 and (fdecl.is_common() or fdecl.is_save()))):
    conv_info.vmap[identifier] = prepend_identifier_if_necessary(identifier)
    ctype = convert_data_type(conv_info=conv_info, fdecl=fdecl, crhs=None)[0]
    rapp = get_rapp()
    rapp("%s& %s = %s;" % (
      ctype, prepend_identifier_if_necessary(identifier), src_var))
  if (crhs is not None):
    return True
  return False

def convert_executable(
      callback, conv_info, args_fdecl_with_dim=None, blockdata=None):
  top_scope = scope(parent=None)
  if (conv_info.fproc.uses_save):
    macro = "FEM_CMN_SVE"
    if (conv_info.fproc.conv_hook.needs_sve_dynamic_parameters):
      macro += "_DYNAMIC_PARAMETERS"
    top_scope.append("%s(%s);" % (macro, conv_info.fproc.name.value))
  top_scope.remember_insert_point()
  curr_scope = top_scope
  if (args_fdecl_with_dim is not None):
    for fdecl in args_fdecl_with_dim:
      declare_size_dim_identifiers(
        conv_info=conv_info,
        top_scope=top_scope,
        curr_scope=top_scope,
        fdecl=fdecl)
      cdims = convert_dims(conv_info=conv_info, dim_tokens=fdecl.dim_tokens)
      top_scope.append("%s(%s);" % (fdecl.id_tok.value, cdims))
  if (blockdata is not None):
    for fproc in blockdata:
      callback("  %s(cmn);" % fproc.name.value)
  if (conv_info.fproc.uses_read):
    top_scope.append("common_read read(cmn);")
  if (conv_info.fproc.uses_write):
    top_scope.append("common_write write(cmn);")
  top_scope_point_before_common = top_scope.current_point()
  for common_name,fdecl_list in conv_info.fproc.common.items():
    if (common_name in conv_info.fproc.conv_hook.variant_common_names):
      continue
    top_scope.remember_insert_point()
    for common_fdecl in fdecl_list:
      fdecl = conv_info.fproc.fdecl_by_identifier.get(
        common_fdecl.id_tok.value)
      if (    fdecl.use_count != 0
          and fdecl.id_tok.value not in conv_info.vmap):
        declare_identifier(
          conv_info=conv_info,
          top_scope=top_scope,
          curr_scope=top_scope,
          id_tok=fdecl.id_tok)
    if (not top_scope.insert_point_is_current()):
      top_scope.top_append("// COMMON %s" % common_name)
  if (not top_scope.point_is_current(point=top_scope_point_before_common)):
    top_scope.append("//")
  top_scope.remember_insert_point()
  def declare_identifiers(id_tokens):
    for id_tok in id_tokens:
      if (id_tok.value not in conv_info.vmap):
        declare_identifier(
          conv_info=conv_info,
          top_scope=top_scope,
          curr_scope=curr_scope,
          id_tok=id_tok)
  from fable.tokenization import extract_identifiers
  variant_buffers = convert_variant_allocate_and_bindings(
    conv_info=conv_info, top_scope=top_scope)
  top_scope.remember_insert_point()
  cei = conv_info.fproc.classified_equivalence_info()
  sve_equivalences = cei.save.equiv_tok_cluster_by_identifier
  for identifier in sorted(conv_info.fproc.fdecl_by_identifier.keys()):
    fdecl = conv_info.fproc.fdecl_by_identifier[identifier]
    if (    fdecl.is_save()
        and fdecl.use_count > 1
        and not identifier in sve_equivalences
        and not identifier in conv_info.vmap):
      declare_identifier(
        conv_info=conv_info,
        top_scope=top_scope,
        curr_scope=top_scope,
        id_tok=fdecl.id_tok)
  if (not top_scope.insert_point_is_current()):
    top_scope.top_append("// SAVE")
    top_scope.append("//")
  if (conv_info.fproc.conv_hook.needs_is_called_first_time):
    first_time_scope = top_scope.open_nested_scope(
      opening_text=["if (is_called_first_time) {"])
    if (len(variant_buffers.first_time) != 0):
      first_time_scope.append(
        "using fem::mbr; // member of variant common or equivalence")
      for mbr_buffer in variant_buffers.first_time:
        mbr_scope = first_time_scope.open_nested_scope(opening_text=["{"])
        for line in mbr_buffer:
          mbr_scope.append(line)
        mbr_scope.close_nested_scope()
    for nlist,clist in conv_info.fproc.data:
      declare_identifiers(
        id_tokens=extract_identifiers(tokens=nlist))
      for repetition_tok,ctoks in clist:
        if (repetition_tok is not None):
          declare_identifiers(
            id_tokens=extract_identifiers(tokens=[repetition_tok]))
        declare_identifiers(
          id_tokens=extract_identifiers(tokens=ctoks))
    if (not conv_info.fproc.conv_hook.data_init_after_variant_bind):
      convert_data(conv_info=conv_info, data_init_scope=first_time_scope)
    first_time_scope.close_nested_scope()
    top_scope.remember_insert_point()
  if (len(variant_buffers.loc_equivalences) != 0):
    top_scope.append("local_equivalences loc_equivalences;")
    mbr_scope = top_scope.open_nested_scope(opening_text=["{"])
    mbr_scope.append("using fem::mbr; // member")
    for line in variant_buffers.loc_equivalences:
      mbr_scope.append(line)
    mbr_scope.close_nested_scope()
  for line in variant_buffers.bindings:
    top_scope.append(line)
  top_scope.remember_insert_point()
  if (conv_info.fproc.conv_hook.data_init_after_variant_bind):
    data_init_scope = top_scope.open_nested_scope(
      opening_text=["if (is_called_first_time) {"])
    convert_data(conv_info=conv_info, data_init_scope=data_init_scope)
    data_init_scope.close_nested_scope()
  top_scope.remember_insert_point()
  from fable.tokenization import fmt_tokens_as_string
  def get_cfmt_from_format(stmt_label):
    fmt_tokens = conv_info.fproc.format.get(stmt_label)
    if (fmt_tokens is None):
      tok.raise_semantic_error(
        "Unknown FORMAT statement label: %s" % tok.value)
    return '"(' + escape_string_literal(fmt_tokens_as_string(
      tokens=fmt_tokens, comma=fmt_comma_placeholder)) + ')"'
  fmt_counts_by_statement_label = \
    conv_info.fproc.fmt_counts_by_statement_label()
  for stmt_label in sorted(fmt_counts_by_statement_label.keys()):
    if (fmt_counts_by_statement_label[stmt_label] > 1):
      cfmt = get_cfmt_from_format(stmt_label=stmt_label)
      top_scope.append(
        "static const char* format_%s = %s;" % (stmt_label, cfmt))
  def curr_scope_append_return_function():
    curr_scope.append(
      "return %s;" % conv_info.vmap[conv_info.fproc.name.value])
  close_scope_after_next_executable = False
  dos_to_close_by_label = {}
  from fable.read import Error
  from fable import SemanticError
  for ei in conv_info.fproc.executable:
    conv_info.comment_manager.insert_before(
      executable_info=ei, callback=curr_scope.append_comment)
    lbl = ei.ssl.label
    if (    lbl is not None
        and lbl in conv_info.fproc.target_statement_labels()
        and not close_scope_after_next_executable):
      curr_scope.append_statement_label(label=lbl)
    def search_for_id_tokens_and_declare_identifiers():
      id_tokens = []
      def callback(tok, next_tok):
        id_tokens.append(tok)
      ei.search_for_id_tokens(callback=callback)
      declare_identifiers(id_tokens=id_tokens)
      return id_tokens
    try:
      if (ei.key == "assignment"):
        lhs_id_tokens = extract_identifiers(tokens=ei.lhs_tokens)
        assert len(lhs_id_tokens) != 0
        rhs_id_tokens = extract_identifiers(tokens=ei.rhs_tokens)
        for id_tokens in lhs_id_tokens[1:], rhs_id_tokens:
          declare_identifiers(id_tokens=id_tokens)
        crhs = convert_tokens(
          conv_info=conv_info, tokens=ei.rhs_tokens)
        id_tok = lhs_id_tokens[0]
        assign_here = id_tok.value in conv_info.vmap
        if (not assign_here):
          assign_here = declare_identifier(
            conv_info=conv_info,
            top_scope=top_scope,
            curr_scope=curr_scope,
            id_tok=id_tok,
            crhs=crhs)
        clhs = convert_tokens(conv_info=conv_info, tokens=ei.lhs_tokens)
        if (assign_here):
          def in_place_op_left():
            if (not crhs.startswith(clhs)): return False
            i = len(clhs)
            if (i == len(crhs)): return False
            if (crhs[i] != " "): return False
            i += 1
            if (i == len(crhs)): return False
            op = crhs[i]
            if (op != "+"): return False
            i += 1
            if (i == len(crhs)): return False
            if (crhs[i] != " "): return False
            i += 1
            if (i == len(crhs)): return False
            if (crhs[i:] == "1"):
              curr_scope.append("%s++;" % clhs)
            else:
              curr_scope.append("%s %s= %s;" % (clhs, op, crhs[i:]))
            return True
          def in_place_op_right():
            if (not crhs.endswith(clhs)): return False
            i = len(crhs) - len(clhs)
            if (i == 0): return False
            i -= 1
            if (crhs[i] != " "): return False
            if (i == 0): return False
            i -= 1
            op = crhs[i]
            if (op != "+"): return False
            if (i == 0): return False
            i -= 1
            if (crhs[i] != " "): return False
            if (i == 0): return False
            if (crhs[:i] == "1"):
              curr_scope.append("%s++;" % clhs)
            else:
              curr_scope.append("%s %s= %s;" % (clhs, op, crhs[:i]))
            return True
          if (not in_place_op_left() and not in_place_op_right()):
            curr_scope.append("%s = %s;" % (clhs, crhs))
      elif (ei.key == "inquire"):
        search_for_id_tokens_and_declare_identifiers()
        iuflist = ei.iuflist
        if (iuflist.unit is not None):
          if (iuflist.file is not None):
            ei.ssl.raise_semantic_error(
              "Conflicting UNIT vs. FILE in INQUIRE statement"
              " (exactly one is needed)", i=ei.start)
          io_function_specialization = "_unit"
          uf_tokens = iuflist.unit
        elif (iuflist.file is not None):
          io_function_specialization = "_file"
          uf_tokens = iuflist.file
        else:
          ei.ssl.raise_semantic_error(
            "Missing UNIT or FILE in INQUIRE statement", i=ei.start)
        io_call_args = convert_tokens(conv_info=conv_info, tokens=uf_tokens)
        convert_io_statement_with_err(
          conv_info=conv_info,
          curr_scope=curr_scope,
          io_function="inquire",
          io_function_specialization=io_function_specialization,
          io_call_args=io_call_args,
          iolist=iuflist)
      elif (ei.key == "file_positioning"):
        search_for_id_tokens_and_declare_identifiers()
        io_call_args = convert_tokens(
          conv_info=conv_info, tokens=ei.alist.unit)
        convert_io_statement_with_err(
          conv_info=conv_info,
          curr_scope=curr_scope,
          io_function=ei.io_function,
          io_function_specialization="",
          io_call_args=io_call_args,
          iolist=ei.alist)
      elif (ei.key == "open"):
        search_for_id_tokens_and_declare_identifiers()
        olist = ei.olist
        if (olist.unit is None):
          ei.ssl.raise_semantic_error(
            "Missing UNIT in OPEN statement", i=ei.start)
        cunit = convert_tokens(conv_info=conv_info, tokens=olist.unit)
        if (olist.file is None):
          cfile = "fem::file_not_specified"
        else:
          cfile = convert_tokens(conv_info=conv_info, tokens=olist.file)
        convert_io_statement_with_err(
          conv_info=conv_info,
          curr_scope=curr_scope,
          io_function="open",
          io_function_specialization="",
          io_call_args="%s, %s" % (cunit, cfile),
          iolist=olist)
      elif (ei.key == "close"):
        search_for_id_tokens_and_declare_identifiers()
        cllist = ei.cllist
        if (cllist.unit is None):
          ei.ssl.raise_semantic_error(
            "Missing UNIT in CLOSE statement", i=ei.start)
        cunit = convert_tokens(conv_info=conv_info, tokens=cllist.unit)
        convert_io_statement_with_err(
          conv_info=conv_info,
          curr_scope=curr_scope,
          io_function="close",
          io_function_specialization="",
          io_call_args=cunit,
          iolist=cllist)
      elif (ei.key in ["read", "write", "print"]):
        search_for_id_tokens_and_declare_identifiers()
        cilist = ei.cilist
        if (ei.key == "print"):
          work_key = "write"
          cunit = "6"
        else:
          work_key = ei.key
          assert cilist.unit is not None
          cunit = convert_tokens(conv_info=conv_info, tokens=cilist.unit)
          if (cunit == "star "): cunit = "6"
        def conv_fmt():
          if (ei.fmt_tokens is not None):
            return '"(' + escape_string_literal(fmt_tokens_as_string(
              tokens=ei.fmt_tokens, comma=fmt_comma_placeholder)) + ')"'
          tl = cilist.fmt
          if (tl is None):
            return None
          if (len(tl) == 1):
            tok = tl[0]
            if (tok.is_op_with(value="*")):
              return "star"
            if (tok.is_integer()):
              stmt_label = tok.value
              if (fmt_counts_by_statement_label[stmt_label] > 1):
                return "format_%s" % stmt_label
              return get_cfmt_from_format(stmt_label=stmt_label)
          return convert_tokens(conv_info=conv_info, tokens=tl)
        cfmt = conv_fmt()
        cchain = []
        has_iostat = False
        for slot in ["rec", "iostat"]:
          tokens = getattr(cilist, slot)
          if (tokens is not None):
            cchain.append("%s(%s)" % (
              slot, convert_tokens(conv_info=conv_info, tokens=tokens)))
            if slot == "iostat":
              has_iostat = True
        if (len(cchain) == 0):
          cchain = ""
        else:
          cchain = "." + ".".join(cchain)
        iolist_id_tokens = extract_identifiers(tokens=ei.iolist)
        declare_identifiers(id_tokens=iolist_id_tokens)
        if (cfmt is None):
          cargs = "%s, fem::unformatted" % cunit
        else:
          cargs = "%s, %s" % (cunit, cfmt)
        implied_dos = []
        find_implied_dos(result=implied_dos, tokens=ei.iolist)
        if (len(implied_dos) == 0):
          if (    cilist.end is None
              and cilist.err is None
              and not has_iostat):
            io_scope = curr_scope
          else:
            io_scope = curr_scope.open_nested_scope(opening_text=["try {"])
          io_op = "%s(%s)%s" % (work_key, cargs, cchain)
          if (len(ei.iolist) == 0):
            io_scope.append(io_op+";")
          else:
            convert_io_loop(
              io_scope=io_scope,
              io_op=io_op,
              conv_info=conv_info,
              tokens=ei.iolist)
        else:
          is_internal_file = False
          if (cilist.unit is not None):
            unit_id_tokens = extract_identifiers(tokens=cilist.unit)
            if (len(unit_id_tokens) >= 1):
              unit_fdecl = conv_info.fproc.get_fdecl(id_tok=unit_id_tokens[0])
              if (    unit_fdecl.data_type is not None
                  and unit_fdecl.data_type.value == "character"):
                is_internal_file = True
          if (   cilist.end is not None
              or cilist.err is not None
              or has_iostat):
            opening_line = "try {"
          else:
            opening_line = "{"
          io_scope = curr_scope.open_nested_scope(opening_text=[opening_line])
          if (is_internal_file): cmn = ""
          else:                  cmn = "cmn, "
          io_scope.append("%s_loop %sloop(%s%s);" % (
            work_key, work_key[0], cmn, cargs))
          if (len(cchain) != 0):
            io_scope.append("%sloop%s;" % (work_key[0], cchain))
          convert_io_loop(
            io_scope=io_scope,
            io_op=work_key[0]+"loop",
            conv_info=conv_info,
            tokens=ei.iolist)
        if (io_scope is not curr_scope):
          io_scope.close_nested_scope()
        for slot,cexception in [("end", "read_end"), ("err", "io_err")]:
          tokens = getattr(cilist, slot)
          if (tokens is not None):
            catch_scope = curr_scope.open_nested_scope(
              opening_text=["catch (fem::%s const&) {" % cexception])
            clabel = convert_tokens(conv_info=conversion_info(), tokens=tokens)
            catch_scope.append("goto statement_%s;" % clabel)
            catch_scope.close_nested_scope()
          else:
            if has_iostat:
              catch_scope = curr_scope.open_nested_scope(
                opening_text=["catch (fem::%s const&) {" % cexception])
              catch_scope.close_nested_scope()
      elif (ei.key == "do"):
        if (ei.id_tok.value not in conv_info.vmap):
          declare_identifier(
            conv_info=conv_info,
            top_scope=top_scope,
            curr_scope=curr_scope,
            id_tok=ei.id_tok)
        for token in ei.tokens:
          declare_identifiers(
            id_tokens=extract_identifiers(tokens=token.value))
        curr_scope = convert_to_fem_do(
          conv_info=conv_info,
          parent_scope=curr_scope,
          i_tok=ei.id_tok,
          fls_tokens=ei.tokens)
        if (ei.label is not None):
          dos_to_close_by_label.setdefault(ei.label, []).append(curr_scope)
      elif (ei.key == "dowhile"):
        declare_identifiers(
          id_tokens=extract_identifiers(tokens=ei.cond_tokens))
        c = convert_tokens(conv_info=conv_info, tokens=ei.cond_tokens)
        curr_scope = curr_scope.open_nested_scope(
          opening_text=["while %s {" % c])
        if (ei.label is not None):
          dos_to_close_by_label.setdefault(ei.label, []).append(curr_scope)
      elif (ei.key == "cycle"):
        curr_scope.append("continue;")
      elif (ei.key == "exit"):
        curr_scope.append("break;")
      elif (ei.key == "enddo"):
        if (dos_to_close_by_label.get(ei.ssl.label) is None):
          curr_scope = curr_scope.close_nested_scope()
      elif (ei.key == "if"):
        declare_identifiers(
          id_tokens=extract_identifiers(tokens=ei.cond_tokens))
        c = convert_tokens(conv_info=conv_info, tokens=ei.cond_tokens)
        curr_scope = curr_scope.open_nested_scope(
          opening_text=["if (%s) {" % c])
        close_scope_after_next_executable = True
        continue
      elif (ei.key == "if_then"):
        declare_identifiers(
          id_tokens=extract_identifiers(tokens=ei.cond_tokens))
        c = convert_tokens(conv_info=conv_info, tokens=ei.cond_tokens)
        curr_scope = curr_scope.open_nested_scope(
          opening_text=["if (%s) {" % c])
      elif (ei.key == "elseif_then"):
        declare_identifiers(
          id_tokens=extract_identifiers(tokens=ei.cond_tokens))
        c = convert_tokens(conv_info=conv_info, tokens=ei.cond_tokens)
        curr_scope = curr_scope.attach_tail(
          opening_text=["else if (%s) {" % c])
      elif (ei.key == "else"):
        curr_scope = curr_scope.attach_tail(opening_text=["else {"])
      elif (ei.key == "endif"):
        curr_scope = curr_scope.close_nested_scope()
      elif (ei.key == "if_arithmetic"):
        declare_identifiers(
          id_tokens=extract_identifiers(tokens=ei.cond_tokens))
        c = convert_tokens(conv_info=conv_info, tokens=ei.cond_tokens)
        curr_scope = curr_scope.open_nested_scope(
          opening_text=["switch (fem::if_arithmetic(%s)) {" % c])
        def lbl(i): return "statement_" + ei.labels[i].value
        curr_scope.append("case -1: goto %s;" % lbl(0))
        curr_scope.append("case  0: goto %s;" % lbl(1))
        curr_scope.append("default: goto %s;" % lbl(2))
        curr_scope = curr_scope.close_nested_scope()
      elif (ei.key == "call"):
        fdecl = conv_info.fproc.get_fdecl(id_tok=ei.subroutine_name)
        if (fdecl.is_intrinsic()):
          from fable import intrinsics
          cmn = ""
          if (ei.subroutine_name.value == "getarg"):
            called = "cmn.getarg"
          elif (ei.subroutine_name.value in intrinsics.io_set_lower):
            called = "cmn.io.%s" % ei.subroutine_name.value
          else:
            called = "fem::%s" % ei.subroutine_name.value
        else:
          if (called_fproc_needs_cmn(
                conv_info=conv_info,
                called_name=ei.subroutine_name.value)):
            cmn = "cmn"
          else:
            cmn = ""
          called = conv_info.vmapped_callable(
            identifier=ei.subroutine_name.value)
        if (ei.arg_token is None):
          curr_scope.append("%s(%s);" % (called, cmn))
        else:
          id_tokens = extract_identifiers(tokens=ei.arg_token.value)
          declare_identifiers(id_tokens=id_tokens)
          a = convert_tokens(
            conv_info=conv_info, tokens=ei.arg_token.value, commas=True)
          def cmn_a():
            if (len(cmn) == 0): return a
            if (len(a) == 0): return cmn
            return cmn + ", " + a
          curr_scope.append("%s(%s);" % (called, cmn_a()))
      elif (ei.key == "return"):
        if (conv_info.fproc.fproc_type == "function"):
          curr_scope_append_return_function()
        elif (ei is not conv_info.fproc.executable[-1]):
          curr_scope.append("return;")
      elif (ei.key == "continue"):
        pass
      elif (ei.key == "goto"):
        curr_scope.append("goto statement_%s;" % ei.label.value)
      elif (ei.key == "goto_computed"):
        search_for_id_tokens_and_declare_identifiers()
        ccond = convert_tokens(conv_info=conv_info, tokens=ei.tokens)
        switch_scope = curr_scope.open_nested_scope(
          opening_text=["switch (%s) {" % ccond])
        for i,label in enumerate(ei.labels):
          switch_scope.append(
            "case %d: goto statement_%s;" % (i+1, label.value))
        switch_scope.append("default: break;")
        switch_scope.close_nested_scope()
      elif (ei.key == "stop"):
        if (ei.arg_token is None):
          cmsg = "0"
        elif (ei.arg_token.is_integer()):
          cmsg = strip_leading_zeros(string=ei.arg_token.value)
        else:
          cmsg = convert_token(vmap={}, leading=True, tok=ei.arg_token)
        curr_scope.append("FEM_STOP(%s);" % cmsg)
      elif (ei.key == "entry"):
        curr_scope.append(
          "// UNHANDLED: ENTRY %s" % ei.ssl.code_with_strings()[5:])
      else:
        curr_scope.append(
          'FEM_THROW_UNHANDLED("executable %s: %s");' % (
            ei.key, ei.ssl.code_with_strings()))
      if (close_scope_after_next_executable):
        close_scope_after_next_executable = False
        curr_scope = curr_scope.close_nested_scope()
      if (ei.ssl.label is not None):
        dos_to_close = dos_to_close_by_label.get(ei.ssl.label)
        if (dos_to_close is not None):
          for do_scope in reversed(dos_to_close):
            curr_scope = do_scope.close_nested_scope()
    except (Error, SemanticError):
      raise
    except Exception:
      print("*"*80)
      print(ei.ssl.format_error(
        i=None,
        msg="Sorry: fable internal error"))
      print("*"*80)
      print()
      raise
  assert curr_scope.parent is None
  if (    conv_info.fproc.fproc_type == "function"
      and len(conv_info.fproc.executable) != 0
      and conv_info.fproc.executable[-1].key != "return"):
    curr_scope_append_return_function()
  conv_info.comment_manager.flush_remaining(
    callback=curr_scope.append_comment)
  curr_scope.finalize()
  curr_scope.collect_text(callback=callback)

def export_save_struct(callback, conv_info):
  cci = conv_info.converted_commons_info
  if (cci is not None):
    buffer = cci.save_struct_buffers.get(conv_info.fproc.name.value)
    if (buffer is not None):
      for line in buffer:
        callback(line)

def produce_fortran_file_comment(conv_info, callback):
  if (conv_info.fortran_file_comments):
    callback("// Fortran file: %s"
      % conv_info.fproc.body_lines[0].source_line_cluster[0].file_name)

def convert_to_cpp_function(
      cpp_callback,
      hpp_callback,
      conv_info,
      declaration_only=False,
      force_not_implemented=False):
  if (not declaration_only):
    export_save_struct(callback=cpp_callback, conv_info=conv_info)
  fptr = []
  cargs = []
  def cargs_append(ctype, name):
    fptr.append(ctype)
    cargs.append(ctype + " " + name)
  if (    conv_info.fproc.needs_cmn
      and not conv_info.fproc.conv_hook.ignore_common_and_save):
    cargs_append("common&", "cmn")
  args_fdecl_with_dim = []
  for id_tok in conv_info.fproc.args:
    if (id_tok.value == "*"):
      cargs_append("fem::star_type const&", "/* UNHANDLED: star argument */")
      continue
    assert id_tok.value not in conv_info.vmap
    fdecl = conv_info.fproc.get_fdecl(id_tok=id_tok)
    conv_info.set_vmap_from_fdecl(fdecl=fdecl)
    assert fdecl.parameter_assignment_tokens is None
    if (fdecl.use_count == 0):
      arg_name = "/* %s */" % prepend_identifier_if_necessary(id_tok.value)
    else:
      arg_name = prepend_identifier_if_necessary(id_tok.value)
    if (    fdecl.data_type is not None
        and fdecl.data_type.value == "character"):
      if (fdecl.dim_tokens is None):
        if cconst(fdecl=fdecl, short=True) == "c":
            cargs_append("const char *", arg_name)
        else:
            cargs_append("str_%sref" % cconst(fdecl=fdecl, short=True), arg_name)
      else:
        if (len(fdecl.dim_tokens) == 1):
          cdim = ""
        else:
          cdim = "%d" % len(fdecl.dim_tokens)
        cargs_append("str_arr_%sref<%s>" % (
          cconst(fdecl=fdecl, short=True), cdim), arg_name)
    elif (not fdecl.is_user_defined_callable()):
      ctype = convert_data_type(conv_info=conv_info, fdecl=fdecl, crhs=None)[0]
      if (fdecl.dim_tokens is None):
#        cargs_append("%s%s&" % (
#          ctype,
#          cconst(fdecl=fdecl, short=False)),
#          prepend_identifier_if_necessary(arg_name))
        cargs_append("%s" % ctype,
          prepend_identifier_if_necessary(arg_name))
      else:
        if (len(fdecl.dim_tokens) == 1):
          t = ctype
        else:
          t = "%s, %d" % (ctype, len(fdecl.dim_tokens))
        if (t.endswith(">")): templs = " "
        else:                 templs = ""
#        cargs_append("arr_%sref<%s%s>" % (
#         cconst(fdecl=fdecl, short=True), t, templs), arg_name)
#        cconst(fdecl=fdecl, short=True), t, templs
        cargs_append("%s *" % ctype, arg_name)
    else:
      passed = conv_info.fproc.externals_passed_by_arg_identifier.get(
        fdecl.id_tok.value)
      if (passed is None or len(passed) == 0):
        ctype = "UNHANDLED"
      else:
        ctype = sorted(passed)[0]
      cargs_append(
        ctype=ctype+"_function_pointer",
        name=arg_name)
    if (fdecl.dim_tokens is not None and fdecl.use_count != 0):
      args_fdecl_with_dim.append(fdecl)
  cdecl = "void"
  if (conv_info.fproc.name is not None):
    fdecl = conv_info.fproc.get_fdecl(id_tok=conv_info.fproc.name)
    if (fdecl.data_type is not None):
      cdecl = convert_data_type(conv_info=conv_info, fdecl=fdecl, crhs=None)[0]
      conv_info.vmap[conv_info.fproc.name.value] = "return_value"
  if (declaration_only):
    cpp_callback("")
    cpp_callback("// forward declaration (dependency cycle)")
    if (conv_info.inline_all):
      cpp_callback("inline")
    cpp_callback("%s %s(%s);" % (
      cdecl,
      prepend_identifier_if_necessary(conv_info.fproc.name.value),
      ", ".join(fptr)))
    return
  if (conv_info.fproc.is_passed_as_external):
    if (hpp_callback is None): cb = cpp_callback
    else:                      cb = hpp_callback
    cb("")
    cb("typedef %s (*%s_function_pointer)(%s);" % (
      cdecl,
      prepend_identifier_if_necessary(conv_info.fproc.name.value),
      ", ".join(fptr)))
  for callback in [hpp_callback, cpp_callback]:
    if (callback is None): continue
    callback("")
    if (callback is cpp_callback):
      produce_leading_comments(callback=callback, fproc=conv_info.fproc)
      produce_fortran_file_comment(conv_info=conv_info, callback=callback)
    if (conv_info.inline_all):
      callback("inline")
    callback(cdecl)
    if (callback is hpp_callback): last = ";"
    else:                          last = ""
    cname = prepend_identifier_if_necessary(conv_info.fproc.name.value)
    if (len(cargs) == 0):
      callback(cname+"()" + last)
    else:
      callback(cname + "(\n  " + ",\n  ".join(cargs) + ")" + last)
  cpp_callback("{")
  if (cdecl != "void"):
    cpp_callback("  %s %s = %s;" % (
      cdecl,
      conv_info.vmap[conv_info.fproc.name.value],
      zero_shortcut_if_possible(ctype=cdecl)))
  if (force_not_implemented):
    cpp_callback("  throw TBXX_NOT_IMPLEMENTED();")
  else:
    convert_executable(
      callback=cpp_callback,
      conv_info=conv_info,
      args_fdecl_with_dim=args_fdecl_with_dim)
  cpp_callback("}")
  produce_trailing_comments(callback=callback, fproc=conv_info.fproc)

def convert_to_struct(
      callback,
      separate_cmn_hpp,
      fproc,
      struct_type,
      struct_name,
      equivalence_simple,
      id_tok_list):
  assert struct_type in ["common", "save"]
  need_dynamic_parameters = False
  conv_info = conversion_info(fproc=fproc)
  callback("")
  callback("struct %s" % struct_name)
  callback("{")
  sve_equivalences = {}
  cmn_equivalences = {}
  have_variant_block = False
  if (    struct_type == "save"
      and conv_info.fproc.conv_hook.needs_is_called_first_time):
    for common_name in sorted(conv_info.fproc.conv_hook.variant_common_names):
      callback("  fem::variant_bindings %s_bindings;" % common_name)
      have_variant_block = True
    #
    cei = conv_info.fproc.classified_equivalence_info()
    sve_equivalences = cei.save.equiv_tok_cluster_by_identifier
    if (len(sve_equivalences) != 0):
      callback("  fem::variant_core_and_bindings save_equivalences;")
      have_variant_block = True
    cmn_equivalences = cei.common.equiv_tok_cluster_by_identifier
  #
  from fable.tokenization import extract_identifiers
  const_identifiers = {}
  const_id_toks = []
  remaining_id_tok_list = []
  for id_tok in id_tok_list:
    if (id_tok.value in sve_equivalences):
      continue
    if (id_tok.value in cmn_equivalences):
      continue
    remaining_id_tok_list.append(id_tok)
    fdecl = conv_info.fproc.get_fdecl(id_tok=id_tok)
    for tokens in [fdecl.size_tokens, fdecl.dim_tokens]:
      if (tokens is None):
        continue
      def parameter_recursion(tokens):
        have_dynamic_dependency = False
        for id_tok in extract_identifiers(tokens=tokens):
          if (id_tok.value in const_identifiers):
            continue
          const_identifiers[id_tok.value] = None
          fdecl = conv_info.fproc.get_fdecl(id_tok=id_tok)
          hdp = parameter_recursion(
            tokens=fdecl.required_parameter_assignment_tokens())
          if (hdp or id_tok.value in conv_info.fproc.dynamic_parameters):
            have_dynamic_dependency = True
          const_identifiers[id_tok.value] = have_dynamic_dependency
          const_id_toks.append(id_tok)
        return have_dynamic_dependency
      parameter_recursion(tokens=tokens)
  initializers = []
  const_definitions = []
    # ISO C++ 9.4.2-4:
    #   "The member shall still be defined in a namespace scope if it is used
    #   in the program and the namespace scope definition shall not contain
    #   an initializer."
  if (len(const_id_toks) != 0):
    if (have_variant_block):
      callback("")
    append_empty_line = False
    for id_tok in const_id_toks:
      fdecl = conv_info.fproc.get_fdecl(id_tok=id_tok)
      ctype = convert_data_type(conv_info=conv_info, fdecl=fdecl, crhs=None)[0]
      if (const_identifiers[id_tok.value]):
        need_dynamic_parameters = True
        callback("  const %s %s;" % (
          ctype, prepend_identifier_if_necessary(id_tok.value)))
        if (id_tok.value in conv_info.fproc.dynamic_parameters):
          crhs = "dynamic_params." + prepend_identifier_if_necessary(
            id_tok.value)
        else:
          crhs = convert_tokens(
            conv_info=conv_info, tokens=fdecl.parameter_assignment_tokens)
        initializers.append((id_tok.value, crhs))
      else:
        crhs = convert_tokens(
          conv_info=conv_info, tokens=fdecl.parameter_assignment_tokens)
        callback("  static const %s %s = %s;" % (
          ctype, prepend_identifier_if_necessary(id_tok.value), crhs))
        const_definitions.append(
          "const %s %s::%s;" % (
            ctype, struct_name, prepend_identifier_if_necessary(id_tok.value)))
        append_empty_line = True
    if (append_empty_line):
      callback("")
  #
  deferred_arr_members = []
  deferred_arr_initializers = []
  for id_tok in remaining_id_tok_list:
    fdecl = conv_info.fproc.get_fdecl(id_tok=id_tok)
    if (fdecl.id_tok.value in const_identifiers):
      continue
    ctype, cdims, crhs = convert_data_type_and_dims(
      conv_info=conv_info, fdecl=fdecl, crhs=None, force_arr=True)[:3]
    if (cdims is None):
      callback("  %s %s;" % (
        ctype, prepend_identifier_if_necessary(id_tok.value)))
      if (crhs is None):
        crhs = zero_shortcut_if_possible(ctype=ctype)
      initializers.append((id_tok.value, crhs))
    elif (not equivalence_simple or need_dynamic_parameters):
      callback("  %s %s;" % (
        ctype, prepend_identifier_if_necessary(id_tok.value)))
      initializers.append((id_tok.value, "%s, fem::fill0" % cdims))
    else:
      ctype_core = convert_data_type(
        conv_info=conv_info, fdecl=fdecl, crhs=None)[0]
      cstatic_size = convert_dims_to_static_size(
        conv_info=conv_info, dim_tokens=fdecl.dim_tokens)
      callback("  %s %s_memory[%s];" % (
        ctype_core,
        prepend_identifier_if_necessary(id_tok.value),
        cstatic_size))
      if (fdecl.data_type.value == "character"):
        deferred_arr_members.append("  str_arr_ref<%d> %s;" % (
          len(fdecl.dim_tokens),
          prepend_identifier_if_necessary(id_tok.value)))
      else:
        deferred_arr_members.append("  %s %s;" % (
          ad_hoc_change_arr_to_arr_ref(ctype=ctype),
          prepend_identifier_if_necessary(id_tok.value)))
      deferred_arr_initializers.append((
        id_tok.value,
        "*%s_memory, %s, fem::fill0" % (
          prepend_identifier_if_necessary(id_tok.value), cdims)))
  if (len(deferred_arr_members) != 0):
    callback("")
    for line in deferred_arr_members:
      callback(line)
    initializers.extend(deferred_arr_initializers)
  n = len(initializers)
  if (n != 0):
    callback("")
    if (not need_dynamic_parameters):
      callback("  %s() :" % struct_name)
    else:
      callback("  %s(" % struct_name)
      callback("    dynamic_parameters const& dynamic_params)")
      callback("  :")
    for i in range(n):
      ii = initializers[i]
      if (i+1 == n): comma = ""
      else:          comma = ","
      callback("    %s(%s)%s" % (
        prepend_identifier_if_necessary(ii[0]), ii[1], comma))
    callback("  {}")
  callback("};")
  #
  if (len(const_definitions) != 0):
    callback("")
    need_ifdef = (separate_cmn_hpp and struct_type == "common")
    if (need_ifdef):
      callback("#ifdef FEM_TRANSLATION_UNIT_WITH_MAIN")
    for cd in const_definitions:
      callback(cd)
    if (need_ifdef):
      callback("#endif")
  #
  return group_args(
    need_dynamic_parameters=need_dynamic_parameters)

def generate_common_report(
      common_fdecl_list_sizes,
      common_equiv_tok_seqs,
      ccode_registry,
      member_registry,
      variant_due_to_equivalence_common_names,
      stringio):
  from six.moves import StringIO
  variant_common_names = set()
  if (stringio is None):
    report = StringIO()
  else:
    report = stringio
  for common_name,fproc_cpp_pairs in ccode_registry.items():
    fprocs_by_cpp = {}
    for fproc,cpp in fproc_cpp_pairs:
      fprocs_by_cpp.setdefault("\n".join(cpp), []).append(fproc)
    if (len(fprocs_by_cpp) != 1):
      variant_common_names.add(common_name)
      fprocs_by_cpp_items = list(fprocs_by_cpp.items())
      def size_key(a):
        return len(a[0])
      fprocs_by_cpp_items.sort(key=size_key, reverse=True)
      import difflib
      diff_function = getattr(difflib, "unified_diff", difflib.ndiff)
      def show_fprocs(label, cpp_fprocs):
        print("procedures %s:" % label,
          " ".join(sorted([fproc.name.value for fproc in cpp_fprocs[1]])), file=report)
      main_cpp_fprocs = fprocs_by_cpp_items[0]
      print("common name:", common_name, file=report)
      print("number of variants:", len(fprocs_by_cpp_items), file=report)
      print("total number of procedures using the common block:",
        sum([len(fprocs) for cpp,fprocs in fprocs_by_cpp_items]), file=report)
      show_fprocs("first", main_cpp_fprocs)
      for other_cpp_fprocs in fprocs_by_cpp_items[1:]:
        show_fprocs("second", other_cpp_fprocs)
        print(" ".join([line for line in diff_function(
                      (main_cpp_fprocs[0]+"\n").splitlines(1),
                      (other_cpp_fprocs[0]+"\n").splitlines(1))]), file=report)
  #
  need_empty_line = False
  for identifier in sorted(member_registry.keys()):
    common_names = member_registry[identifier]
    if (len(common_names) != 1):
      print("Name clash: %s in COMMONs: %s" % (
        identifier, ", ".join(sorted(common_names))), file=report)
      need_empty_line = True
  if (need_empty_line):
    print(file=report)
  #
  vv = list(variant_due_to_equivalence_common_names - variant_common_names)
  if (len(vv) != 0):
    print("common variants due to equivalence:", len(vv), file=report)
    size_sums = {}
    for common_name,sizes in common_fdecl_list_sizes.items():
      size_sums[common_name] = sum(sizes)

    vv.sort(key=lambda element: (-size_sums[element], element))
    print("  %-20s   procedures    sum of members" % "common name", file=report)
    for common_name in vv:
      print("  %-20s   %8d         %8d" % (
        common_name,
        len(common_fdecl_list_sizes[common_name]),
        size_sums[common_name]), file=report)
    print(file=report)
    print("Locations of equivalence statements:", file=report)
    reported_already = set()
    for common_name in vv:
      print("  %s" % common_name, file=report)
      prev_loc = ""
      tab = []
      max_len_col1 = 6
      for tok_seq in common_equiv_tok_seqs[common_name]:
        sl, i = tok_seq.stmt_location()
        tag = (sl.file_name, sl.line_number, i)
        if (tag in reported_already):
          break
        reported_already.add(tag)
        vn = tok_seq.value[0].value
        dn, bn = os.path.split(sl.file_name)
        loc = ("%s(%s) %s" % (bn, sl.line_number, dn)).rstrip()
        if (loc == prev_loc): loc = ""
        else: prev_loc = loc
        tab.append((vn, loc))
        max_len_col1 = max(max_len_col1, len(vn))
      if (len(tab) != 0):
        fmt = "    %%-%ds %%s" % max_len_col1
        for row in tab:
          print(fmt % row, file=report)
  #
  if (len(report.getvalue()) != 0 and stringio is None):
    import sys
    report_file_name = "fable_cout_common_report"
    from libtbx.str_utils import show_string
    print("Writing file:", show_string(report_file_name), file=sys.stderr)
    open(report_file_name, "w").write(report.getvalue())
  #
  return variant_common_names

def convert_commons(
      callback,
      separate_cmn_hpp,
      topological_fprocs,
      dynamic_parameters,
      common_equivalence_simple,
      common_report_stringio):
  if (dynamic_parameters is not None):
    callback("")
    callback("struct dynamic_parameters")
    callback("{")
    for dp_props in dynamic_parameters:
      callback("  %s %s;" % (dp_props.ctype, dp_props.name))
    callback("""
  dynamic_parameters(
    fem::command_line_arguments const& command_line_args)
  :""")
    for dp_props in dynamic_parameters:
      if (dp_props is not dynamic_parameters[-1]): c = ","
      else:                                        c = ""
      callback("    %s(%s)%s" % (
        prepend_identifier_if_necessary(dp_props.name),
        str(dp_props.default),
        c))
    callback("""\
  {
    fem::dynamic_parameters_from(command_line_args, %d)"""
      % len(dynamic_parameters))
    for dp_props in dynamic_parameters:
      callback("      .reset_if_given(%s)"
        % prepend_identifier_if_necessary(dp_props.name))
    callback("    ;")
    callback("  }")
    callback("};")
    callback("")
    callback("typedef")
    callback("  fem::dynamic_parameters_capsule<dynamic_parameters>")
    callback("    dynamic_parameters_capsule;")
  #
  common_fdecl_list_sizes = {}
  common_equiv_tok_seqs = {}
  common_ccode_registry = {}
  member_registry = {}
  variant_common_names = set()
  bottom_up_filtered = []
  for fproc in topological_fprocs.bottom_up_list:
    if (not fproc.conv_hook.ignore_common_and_save):
      bottom_up_filtered.append(fproc)
  struct_commons_need_dynamic_parameters = set()
  for fproc in bottom_up_filtered:
    fproc.conv_hook.needs_variant_bind = False
    for common_name,common_fdecl_list in fproc.common.items():
      common_fdecl_list_sizes.setdefault(common_name, []).append(
        len(common_fdecl_list))
      id_tok_list = []
      for common_fdecl in common_fdecl_list:
        assert common_fdecl.size_tokens is None
        id_tok_list.append(common_fdecl.id_tok)
        member_registry.setdefault(
          common_fdecl.id_tok.value, set()).add(common_name)
        if (common_name not in common_equivalence_simple):
          equiv_tok_cluster = fproc.equivalence_info() \
            .equiv_tok_cluster_by_identifier.get(common_fdecl.id_tok.value)
          if (equiv_tok_cluster is not None):
            fproc.conv_hook.needs_variant_bind = True
            variant_common_names.add(common_name)
            for equiv_tok in equiv_tok_cluster:
              for tok_seq in equiv_tok.value:
                common_equiv_tok_seqs.setdefault(common_name, []).append(
                  tok_seq)
      struct_name = "common_" + common_name
      buffer = []
      info = convert_to_struct(
        callback=buffer.append,
        separate_cmn_hpp=separate_cmn_hpp,
        fproc=fproc,
        struct_type="common",
        struct_name=struct_name,
        equivalence_simple=(common_name in common_equivalence_simple),
        id_tok_list=id_tok_list)
      if (info.need_dynamic_parameters):
        struct_commons_need_dynamic_parameters.add(struct_name)
      common_ccode_registry.setdefault(common_name, []).append(
        (fproc, buffer))
  variant_common_names.update(generate_common_report(
    common_fdecl_list_sizes=common_fdecl_list_sizes,
    common_equiv_tok_seqs=common_equiv_tok_seqs,
    ccode_registry=common_ccode_registry,
    member_registry=member_registry,
    variant_due_to_equivalence_common_names=variant_common_names,
    stringio=common_report_stringio))
  commons_defined_already = set()
  struct_commons = []
  variant_commons = []
  for fproc in bottom_up_filtered:
    fproc.conv_hook.variant_common_names = set()
    for common_name,common_fdecl_list in fproc.common.items():
      if (common_name in variant_common_names):
        fproc.conv_hook.variant_common_names.add(common_name)
        if (common_name not in commons_defined_already):
          commons_defined_already.add(common_name)
          variant_commons.append(common_name)
      else:
        if (common_name not in commons_defined_already):
          commons_defined_already.add(common_name)
          struct_commons.append("common_"+common_name)
          for line in common_ccode_registry[common_name][0][1]:
            callback(line)
  #
  for fproc in bottom_up_filtered:
    if (not fproc.conv_hook.needs_variant_bind):
      fproc.conv_hook.needs_variant_bind = (
           len(fproc.conv_hook.variant_common_names) != 0
        or fproc.classified_equivalence_info().has_save())
    fproc.conv_hook.needs_is_called_first_time = (
         fproc.conv_hook.needs_variant_bind
      or len(fproc.data) != 0)
    fproc.conv_hook.data_init_after_variant_bind = (
          fproc.conv_hook.needs_variant_bind
      and len(fproc.data) != 0)
    if (fproc.conv_hook.needs_is_called_first_time):
      fproc.uses_save = True
  topological_fprocs.each_fproc_update_is_modified()
  topological_fprocs.each_fproc_update_needs_cmn()
  #
  save_struct_buffers = {}
  save_struct_names = []
  for fproc in bottom_up_filtered:
    id_tok_list = []
    for fdecl in fproc.fdecl_by_identifier.values():
      if (fdecl.is_save()):
        id_tok_list.append(fdecl.id_tok)
    if (    len(id_tok_list) == 0
        and not fproc.conv_hook.needs_is_called_first_time):
      continue
    id_tok_list.sort(key=lambda token: token.value)
    struct_name = "%s_save" % fproc.name.value
    buffer = []
    info = convert_to_struct(
      callback=buffer.append,
      separate_cmn_hpp=separate_cmn_hpp,
      fproc=fproc,
      struct_type="save",
      struct_name=struct_name,
      equivalence_simple=False,
      id_tok_list=id_tok_list)
    save_struct_buffers[fproc.name.value] = buffer
    if (info.need_dynamic_parameters):
      fproc.conv_hook.needs_sve_dynamic_parameters = True
    save_struct_names.append(struct_name)
  if (    len(commons_defined_already) == 0
      and len(save_struct_names) == 0
      and dynamic_parameters is None):
    callback("")
    callback("using fem::common;")
    return
  callback("")
  callback("struct common :")
  leading_bases = ["fem::common"]
  if (dynamic_parameters is not None):
    leading_bases.append("dynamic_parameters_capsule")
  callback("  " + ",\n  ".join(leading_bases + struct_commons))
  callback("{")
  need_empty_line = False
  for common_name in variant_commons:
    callback("  fem::variant_core common_%s;" % common_name)
    need_empty_line = True
  def save_as_sve(struct_name): return struct_name[:-3]+"ve"
  for struct_name in save_struct_names:
    callback("  fem::cmn_sve %s;" % save_as_sve(struct_name))
    need_empty_line = True
  if (need_empty_line):
    callback("")
  initializations = ["fem::common(argc, argv)"]
  if (dynamic_parameters is not None):
    initializations.append(
      "dynamic_parameters_capsule(command_line_args)")
    for struct_name in struct_commons:
      if (struct_name in struct_commons_need_dynamic_parameters):
        initializations.append("%s(dynamic_params)" % struct_name)
  callback("""\
  common(
    int argc,
    char const* argv[])
  :
    %s
  {}""" % ",\n    ".join(initializations))
  callback("};")
  #
  return group_args(
    member_registry=member_registry,
    save_struct_buffers=save_struct_buffers)

include_fem_hpp = \
  "#include <fem.hpp> // Fortran EMulation library of fable module"

def include_guard(callback, namespace, suffix):
  s = namespace.upper().replace("::", "_") + suffix
  callback("#ifndef %s" % s)
  callback("#define %s" % s)
  callback("")

def open_namespace(callback, namespace, using_namespace_major_types=True):
  ns = namespace.split("::")
  for component in ns:
    callback("namespace %s {" % component)
  if (using_namespace_major_types):
    callback("""
using namespace fem::major_types;""")
  return ns

def close_namespace(callback, namespace, hpp_guard):
  callback("")
  ns = namespace.split("::")
  callback("%s // namespace %s" % ("}"*len(ns), namespace))
  if (hpp_guard):
    callback("")
    callback("#endif // GUARD")
  return ns

class hpp_cpp_buffers(object):

  __slots__ = ["hpp", "cpp"]

  def __init__(O):
    O.hpp = []
    O.cpp = []

def convert_program(callback, global_conv_info, namespace, hpp_guard, debug):
  main_calls = []
  for fproc in global_conv_info.topological_fprocs.bottom_up_list:
    if (not fproc.is_program()): continue
    conv_info = global_conv_info.specialized(fproc=fproc)
    export_save_struct(callback=callback, conv_info=conv_info)
    cname = fproc.name.value
    main_calls.append(cname)
    callback("")
    produce_leading_comments(callback=callback, fproc=fproc)
    produce_fortran_file_comment(conv_info=conv_info, callback=callback)
    callback("""\
void
%s(
  int argc,
  char const* argv[])
{""" % cname)
    if (not fproc.needs_cmn):
      callback("""\
  if (argc != 1) {
    throw std::runtime_error("Unexpected command-line arguments.");
  }""")
    result_buffer = []
    try:
      convert_executable(
        callback=result_buffer.append,
        conv_info=conv_info,
        blockdata=global_conv_info.topological_fprocs.all_fprocs.blockdata)
    except Exception:
      if (not debug): raise
      show_traceback()
    else:
      if (fproc.needs_cmn and not fproc.conv_hook.ignore_common_and_save):
        callback("  common cmn(argc, argv);")
      for line in result_buffer:
        callback(line)
      callback("}")
    produce_trailing_comments(callback=callback, fproc=fproc)
  #
  ns = close_namespace(
    callback=callback, namespace=namespace, hpp_guard=hpp_guard)
  #
  if (len(main_calls) != 0):
    callback("")
    callback("""\
int
main(
  int argc,
  char const* argv[])
{
  return fem::main_with_catch(
    argc, argv,
    %s);
}""" % "::".join(ns + [main_calls[0]]))

def get_missing_external_return_type(fdecls):
  for fdecl in fdecls:
    if (fdecl.data_type is not None):
      return convert_data_type(
        conv_info=conversion_info(), fdecl=fdecl, crhs=None)[0]
  return "void"

default_arr_nd_size_max = 256

def process(
      file_names=None,
      all_fprocs=None,
      top_procedures=None,
      namespace="please_specify",
      include_prefix=None,
      include_guard_suffix=None,
      top_cpp_file_name=None,
      dynamic_parameters=None,
      fortran_file_comments=False,
      fem_do_safe=True,
      arr_nd_size_max=default_arr_nd_size_max,
      inline_all=False,
      common_equivalence_simple=set(),
      suppress_program=False,
      suppress_common=False,
      separate_cmn_hpp=False,
      number_of_function_files=None,
      separate_files_main_namespace={},
      write_separate_files_main_namespace="All",
      separate_files_separate_namespace={},
      write_separate_files_separate_namespace="All",
      ignore_common_and_save=set(),
      force_not_implemented=set(),
      ignore_missing=set(),
      suppress_functions=set(),
      suppress_function_definitions=set(),
      common_report_stringio=None,
      data_values_block_size=8,
      data_specializations=True,
      debug=False):
  assert [file_names, all_fprocs].count(None) == 1
  if (namespace is None or namespace == "please_specify"):
    namespace = "placeholder_please_replace"
  import fable.read
  if (all_fprocs is None):
    all_fprocs = fable.read.process(file_names=file_names)
  for fproc in all_fprocs.all_in_input_order:
    fproc.conv_hook = conv_hook_info()
    fproc.conv_hook.ignore_common_and_save = (
      fproc.name.value in ignore_common_and_save)
  result = []
  def callback(line):
    if (len(result) == 0): prev_line = None
    else:                  prev_line = result[-1]
    lines = break_lines(cpp_text=[line+"\n"], prev_line=prev_line)
    if (len(lines) != 0):
      if (debug):
        print("\n".join(lines))
      result.extend(lines)
  #
  need_function_hpp = False
  if (len(separate_files_main_namespace) != 0):
    need_function_hpp = True
  if (number_of_function_files is not None):
    assert number_of_function_files > 0
    need_function_hpp = True
  if (need_function_hpp):
    separate_cmn_hpp = True
  #
  if (include_guard_suffix is not None):
    include_guard(
      callback=callback, namespace=namespace, suffix=include_guard_suffix)
  #
  if (separate_cmn_hpp):
    callback("#define FEM_TRANSLATION_UNIT_WITH_MAIN")
    callback("")
  #
  def include_separate(callback):
    if (len(separate_files_separate_namespace) == 0):
      return False
    for name in sorted(separate_files_separate_namespace.keys()):
      callback('#include "%s.hpp"' % name)
    return True
  #
  def include_with_prefix(name):
    if (include_prefix is None):
      return '#include "%s.hpp"' % name
    return '#include <%s/%s.hpp>' % (include_prefix, name)
  #
  need_using_major_types = False
  if (need_function_hpp):
    callback(include_with_prefix("functions"))
  elif (separate_cmn_hpp):
    callback(include_with_prefix("cmn"))
  else:
    callback(include_fem_hpp)
    need_using_major_types = True
  callback("")
  if (not need_function_hpp):
    if (include_separate(callback=callback)):
      callback("")
  open_namespace(
    callback=callback,
    namespace=namespace,
    using_namespace_major_types=need_using_major_types)
  #
  topological_fprocs = all_fprocs.build_bottom_up_fproc_list_following_calls(
    top_procedures=top_procedures)
  missing = topological_fprocs.missing_external_fdecls_by_identifier
  if (len(missing) != 0):
    for identifier in sorted(missing.keys()):
      if (identifier in ignore_missing):
        continue
      return_type = get_missing_external_return_type(
        fdecls=missing[identifier])
      callback("""
%s
%s(...)
{
  throw std::runtime_error(
    "Missing function implementation: %s");
}""" % (return_type, identifier, identifier))
  #
  dep_cycles = topological_fprocs.dependency_cycles
  if (len(dep_cycles) != 0):
    callback("")
    callback("/* Dependency cycles: %d" % len(dep_cycles))
    for cycle in dep_cycles:
      callback("     " + " ".join(cycle))
    callback(" */")
  #
  if (dynamic_parameters is not None):
    assert len(dynamic_parameters) != 0
    for fproc in topological_fprocs.bottom_up_list:
      for dp_props in dynamic_parameters:
        fdecl = fproc.fdecl_by_identifier.get(dp_props.name)
        if (fdecl is not None):
          fproc.dynamic_parameters.add(dp_props.name)
  #
  if (separate_cmn_hpp):
    cmn_buffer = []
    cmn_callback = cmn_buffer.append
    include_guard(
      callback=cmn_callback, namespace=namespace, suffix="_CMN_HPP")
    cmn_callback(include_fem_hpp)
    cmn_callback("")
    open_namespace(callback=cmn_callback, namespace=namespace)
  elif (suppress_common):
    def cmn_callback(line): pass
  else:
    cmn_callback = callback
  try:
    converted_commons_info = convert_commons(
      callback=cmn_callback,
      separate_cmn_hpp=separate_cmn_hpp,
      topological_fprocs=topological_fprocs,
      dynamic_parameters=dynamic_parameters,
      common_equivalence_simple=common_equivalence_simple,
      common_report_stringio=common_report_stringio)
  except Exception:
    if (not debug): raise
    show_traceback()
    common_commons_info = None
  if (separate_cmn_hpp):
    close_namespace(callback=cmn_callback, namespace=namespace, hpp_guard=True)
    with open("cmn.hpp", "w") as f:
      print("\n".join(break_lines(cpp_text=cmn_buffer)), file=f)
  #
  separate_function_buffers = []
  separate_function_buffer_by_function_name = {}
  for name,identifiers in separate_files_main_namespace.items():
    if (len(identifiers) == 0):
      raise RuntimeError(
        "separate_files_main_namespace: empty list: %s" % name)
    buffer = []
    buffer.append(include_with_prefix("functions"))
    buffer.append("")
    separate_function_buffers.append((name, buffer))
    open_namespace(callback=buffer.append, namespace=namespace)
    for identifier in identifiers:
      if (identifier in separate_function_buffer_by_function_name):
        raise RuntimeError(
          "separate_files_main_namespace:"
          " ambiguous assignment: %s" % identifier)
      separate_function_buffer_by_function_name[identifier] = buffer
  #
  separate_namespaces = {}
  separate_namespaces_buffers = {}
  for name,identifiers in separate_files_separate_namespace.items():
    if (len(identifiers) == 0):
      raise RuntimeError(
        "separate_files_separate_namespace: empty list: %s" % name)
    buffers = hpp_cpp_buffers()
    for ext in ["hpp", "cpp"]:
      buffer = getattr(buffers, ext)
      if (ext == "hpp"):
        include_guard(callback=buffer.append, namespace=name, suffix="_HPP")
        buffer.append(include_fem_hpp)
      else:
        buffer.append('#include "%s.hpp"' % name)
      buffer.append("")
      open_namespace(callback=buffer.append, namespace=name)
    for identifier in identifiers:
      if (identifier in separate_namespaces):
        raise RuntimeError(
          "separate_files_separate_namespace:"
          " ambiguous assignment: %s" % identifier)
      separate_namespaces[identifier] = name
      separate_namespaces_buffers[identifier] = buffers
  #
  if (not need_function_hpp):
    function_declarations = None
    function_definitions = None
  else:
    function_declarations = []
    function_definitions = []
  #
  global_conv_info = global_conversion_info(
    topological_fprocs=topological_fprocs,
    dynamic_parameters=dynamic_parameters,
    fortran_file_comments=fortran_file_comments,
    fem_do_safe=fem_do_safe,
    arr_nd_size_max=arr_nd_size_max,
    inline_all=inline_all,
    converted_commons_info=converted_commons_info,
    separate_namespaces=separate_namespaces,
    data_values_block_size=data_values_block_size,
    data_specializations=data_specializations)
  #
  for fproc in topological_fprocs.bottom_up_list:
    if (fproc.is_program()):
      continue
    if (fproc.name.value in suppress_functions):
      continue
    hpp_callback = None
    cpp_callback = None
    suppress_cpp = (fproc.name.value in suppress_function_definitions)
    buffers = separate_namespaces_buffers.get(fproc.name.value)
    if (buffers is None):
      if (not need_function_hpp):
        if (not suppress_cpp):
          cpp_callback = callback
      else:
        function_hpp_buffer = []
        function_declarations.append(function_hpp_buffer)
        hpp_callback = function_hpp_buffer.append
        if (not suppress_cpp):
          buffer = separate_function_buffer_by_function_name.get(
            fproc.name.value)
          if (buffer is None):
            if (number_of_function_files is None):
              cpp_callback = callback
            else:
              function_cpp_buffer = []
              function_definitions.append(function_cpp_buffer)
              cpp_callback = function_cpp_buffer.append
          else:
            cpp_callback = buffer.append
    else:
      hpp_callback = buffers.hpp.append
      if (not suppress_cpp):
        cpp_callback = buffers.cpp.append
    if (cpp_callback is None):
      cpp_diverted = []
      cpp_callback = cpp_diverted.append
      if (hpp_callback is None):
        hpp_callback = callback
    if (not need_function_hpp):
      fwds = topological_fprocs.forward_uses_by_identifier.get(
        fproc.name.value)
      if (fwds is not None):
        for fwd_identifier in fwds:
          fwd_fproc = all_fprocs.fprocs_by_name()[fwd_identifier]
          try:
            convert_to_cpp_function(
              hpp_callback=None,
              cpp_callback=cpp_callback,
              conv_info=global_conv_info.specialized(fproc=fwd_fproc),
              declaration_only=True)
          except Exception:
            if (not debug): raise
            show_traceback()
    try:
      convert_to_cpp_function(
        hpp_callback=hpp_callback,
        cpp_callback=cpp_callback,
        conv_info=global_conv_info.specialized(fproc=fproc),
        force_not_implemented=(fproc.name.value in force_not_implemented))
    except Exception:
      if (not debug): raise
      show_traceback()
  #
  for name,buffer in separate_function_buffers:
    close_namespace(
      callback=buffer.append, namespace=namespace, hpp_guard=False)
    if (write_separate_files_main_namespace == "All"
          or name in write_separate_files_main_namespace):
      with open(name+".cpp", "w") as f:
        print("\n".join(break_lines(cpp_text=buffer)), file=f)
  #
  for name,identifiers in separate_files_separate_namespace.items():
    buffers = separate_namespaces_buffers[identifiers[0]]
    for ext in ["hpp", "cpp"]:
      buffer = getattr(buffers, ext)
      close_namespace(
        callback=buffer.append, namespace=name, hpp_guard=(ext=="hpp"))
      if (write_separate_files_separate_namespace == "All"
            or name in write_separate_files_separate_namespace):
        with open(name+"."+ext, "w") as f:
          print("\n".join(break_lines(cpp_text=buffer)), file=f)
  #
  if (function_declarations is not None):
    def write_functions(buffers, serial=None):
      if (buffers is function_declarations):
        assert serial is None
        fn = "functions.hpp"
      elif (serial is None):
        fn = "functions.cpp"
      else:
        fn = "functions_%03d.cpp" % serial
      f = open(fn, "w")
      def fcb(line): print(line, file=f)
      if (buffers is function_declarations):
        include_guard(
          callback=fcb, namespace=namespace, suffix="_FUNCTIONS_HPP")
        fcb(include_with_prefix("cmn"))
        include_separate(callback=fcb)
      else:
        fcb(include_with_prefix("functions"))
      fcb("")
      open_namespace(
        callback=fcb, namespace=namespace, using_namespace_major_types=False)
      for lines in buffers:
        for line in break_lines(cpp_text=lines): fcb(line)
      close_namespace(
        callback=fcb,
        namespace=namespace,
        hpp_guard=(buffers is function_declarations))
      f.close()
    write_functions(function_declarations)
    if (function_definitions is not None and len(function_definitions) != 0):
      buffer_blocks = create_buffer_blocks(
        target_number_of_blocks=number_of_function_files,
        buffers=function_definitions)
      if (len(buffer_blocks) == 1):
        write_functions(buffers=buffer_blocks[0])
      else:
        serial = 0
        for buffers in buffer_blocks:
          serial += 1
          write_functions(buffers=buffers, serial=serial)
  #
  hpp_guard = (include_guard_suffix is not None)
  if (suppress_program):
    close_namespace(
      callback=callback, namespace=namespace, hpp_guard=hpp_guard)
  else:
    try:
      convert_program(
        callback=callback,
        global_conv_info=global_conv_info,
        namespace=namespace,
        hpp_guard=hpp_guard,
        debug=debug)
    except Exception:
      if (not debug): raise
      show_traceback()
  #
  if (top_cpp_file_name is not None):
    with open(top_cpp_file_name, "w") as f:
      print("\n".join(result), file=f)
  #
  return result
