# Licensed under a 3-clause BSD style license - see LICENSE.rst


# This file is automatically generated. Do not edit.
_tabversion = '3.8'

_lr_method = 'LALR'

_lr_signature = 'A63D4C561E2ED1A045DB279536CAFDDA'
    
_lr_action_items = {'CARET':([16,17,29,42,],[39,39,39,39,]),'FUNCNAME':([0,1,4,5,7,8,9,14,15,16,17,19,20,21,22,24,25,27,31,33,34,35,36,45,52,55,58,59,65,66,67,69,72,73,74,76,77,79,81,],[18,-16,18,-32,-17,18,18,-14,-48,-37,-20,-15,18,-33,18,18,-46,-47,18,18,-53,-52,-36,-21,-18,-34,-38,-35,-23,-27,-22,-26,-51,-19,-55,-39,-24,-28,-25,]),'CLOSE_PAREN':([1,2,4,5,7,9,10,12,14,16,17,19,21,23,26,28,30,32,34,35,36,45,47,48,49,50,51,52,55,56,57,58,59,61,62,63,65,66,67,69,71,72,73,74,75,76,77,78,79,81,83,],[-16,-4,-10,-32,-17,-31,-1,-7,-14,-37,-20,-15,-33,-5,-8,-2,55,-30,-53,-52,-36,-21,-13,-11,-6,-9,-3,-18,-34,-29,74,-38,-35,-42,-41,76,-23,-27,-22,-26,79,-51,-19,-55,-40,-39,-24,81,-28,-25,-43,]),'UFLOAT':([0,6,13,33,41,60,68,70,],[-45,-44,34,-45,-45,34,-44,-45,]),'$end':([1,2,3,4,5,7,9,10,12,14,16,17,19,21,23,26,28,32,34,35,36,45,47,48,49,50,51,52,55,56,58,59,65,66,67,69,72,73,74,76,77,79,81,],[-16,-4,0,-10,-32,-17,-31,-1,-7,-14,-37,-20,-15,-33,-5,-8,-2,-30,-53,-52,-36,-21,-13,-11,-6,-9,-3,-18,-34,-29,-38,-35,-23,-27,-22,-26,-51,-19,-55,-39,-24,-28,-25,]),'SOLIDUS':([0,1,2,4,5,7,9,10,14,16,17,19,21,23,24,25,27,28,32,33,34,35,36,45,48,49,51,52,55,56,58,59,65,66,67,69,72,73,74,75,76,77,79,81,],[15,-16,15,15,-32,-17,-31,-12,-14,-37,-20,-15,-33,15,15,-46,-47,-12,-30,15,-53,-52,-36,-21,-11,15,-12,-18,-34,-29,-38,-35,-23,-27,-22,-26,-51,-19,-55,15,-39,-24,-28,-25,]),'UINT':([0,6,7,13,15,16,17,33,34,35,37,38,39,40,41,43,44,53,54,60,64,68,70,80,82,],[17,-44,29,35,-48,-45,42,17,-53,-52,58,-45,-50,-49,-45,66,-45,72,-45,75,-45,72,-45,-45,83,]),'SIGN':([0,15,16,17,29,33,38,39,40,41,42,44,46,54,64,70,80,],[6,-48,6,43,53,6,6,-50,-49,6,53,68,53,6,6,68,6,]),'DOUBLE_STAR':([16,17,29,42,],[40,40,40,40,]),'OPEN_PAREN':([0,1,4,5,7,8,9,11,14,15,16,17,18,19,20,21,22,24,25,27,31,33,34,35,36,38,39,40,44,45,52,54,55,58,59,64,65,66,67,69,72,73,74,76,77,79,81,],[8,-16,8,-32,-17,8,8,33,-14,-48,41,46,-54,-15,8,-33,8,8,-46,-47,8,8,-53,-52,-36,41,-50,-49,70,-21,-18,41,-34,-38,-35,41,-23,-27,-22,-26,-51,-19,-55,-39,-24,-28,-25,]),'PERIOD':([1,4,5,7,9,14,16,17,19,21,34,35,36,45,52,55,58,59,65,66,67,69,72,73,74,76,77,79,81,],[-16,27,-32,-17,27,-14,-37,-20,-15,-33,-53,-52,-36,-21,-18,-34,-38,-35,-23,-27,-22,-26,-51,-19,-55,-39,-24,-28,-25,]),'STAR':([1,4,5,7,9,14,16,17,19,21,34,35,36,45,52,55,58,59,65,66,67,69,72,73,74,76,77,79,81,],[-16,25,-32,-17,25,-14,-37,-20,-15,-33,-53,-52,-36,-21,-18,-34,-38,-35,-23,-27,-22,-26,-51,-19,-55,-39,-24,-28,-25,]),'UNIT':([0,1,4,5,7,8,9,14,15,16,17,19,20,21,22,24,25,27,31,33,34,35,36,45,52,55,58,59,65,66,67,69,72,73,74,76,77,79,81,],[16,-16,16,-32,-17,16,16,-14,-48,-37,-20,-15,16,-33,16,16,-46,-47,16,16,-53,-52,-36,-21,-18,-34,-38,-35,-23,-27,-22,-26,-51,-19,-55,-39,-24,-28,-25,]),}

_lr_action = {}
for _k, _v in _lr_action_items.items():
   for _x,_y in zip(_v[0],_v[1]):
      if not _x in _lr_action:  _lr_action[_x] = {}
      _lr_action[_x][_k] = _y
del _lr_action_items

_lr_goto_items = {'power':([16,17,29,42,],[38,44,54,64,]),'sign':([0,16,33,38,41,44,54,64,70,80,],[13,37,13,37,60,37,37,37,60,82,]),'factor_int':([0,33,],[1,1,]),'numeric_power':([16,38,44,54,64,],[36,59,67,73,77,]),'division_product_of_units':([0,4,24,33,],[2,23,49,2,]),'signed_int':([17,29,42,44,46,70,],[45,52,65,69,71,78,]),'signed_float':([0,33,41,70,],[7,7,62,62,]),'factor':([0,33,],[4,4,]),'function':([0,4,8,9,20,22,24,31,33,],[5,5,5,5,5,5,5,5,5,]),'factor_fits':([0,33,],[14,14,]),'product':([4,9,],[24,31,]),'paren_expr':([41,70,],[63,63,]),'main':([0,33,],[3,57,]),'inverse_unit':([0,4,24,33,],[12,26,50,12,]),'factor_float':([0,33,],[19,19,]),'division':([0,2,4,23,24,33,49,75,],[20,22,20,22,20,20,22,80,]),'unit_expression':([0,4,8,9,20,22,24,31,33,],[9,9,9,9,47,9,9,9,9,]),'product_of_units':([0,4,8,9,22,24,31,33,],[10,28,30,32,48,51,56,10,]),'function_name':([0,4,8,9,20,22,24,31,33,],[11,11,11,11,11,11,11,11,11,]),'frac':([41,70,],[61,61,]),'unit_with_power':([0,4,8,9,20,22,24,31,33,],[21,21,21,21,21,21,21,21,21,]),}

_lr_goto = {}
for _k, _v in _lr_goto_items.items():
   for _x, _y in zip(_v[0], _v[1]):
       if not _x in _lr_goto: _lr_goto[_x] = {}
       _lr_goto[_x][_k] = _y
del _lr_goto_items
_lr_productions = [
  ("S' -> main","S'",1,None,None,None),
  ('main -> product_of_units','main',1,'p_main','generic.py',186),
  ('main -> factor product_of_units','main',2,'p_main','generic.py',187),
  ('main -> factor product product_of_units','main',3,'p_main','generic.py',188),
  ('main -> division_product_of_units','main',1,'p_main','generic.py',189),
  ('main -> factor division_product_of_units','main',2,'p_main','generic.py',190),
  ('main -> factor product division_product_of_units','main',3,'p_main','generic.py',191),
  ('main -> inverse_unit','main',1,'p_main','generic.py',192),
  ('main -> factor inverse_unit','main',2,'p_main','generic.py',193),
  ('main -> factor product inverse_unit','main',3,'p_main','generic.py',194),
  ('main -> factor','main',1,'p_main','generic.py',195),
  ('division_product_of_units -> division_product_of_units division product_of_units','division_product_of_units',3,'p_division_product_of_units','generic.py',207),
  ('division_product_of_units -> product_of_units','division_product_of_units',1,'p_division_product_of_units','generic.py',208),
  ('inverse_unit -> division unit_expression','inverse_unit',2,'p_inverse_unit','generic.py',218),
  ('factor -> factor_fits','factor',1,'p_factor','generic.py',224),
  ('factor -> factor_float','factor',1,'p_factor','generic.py',225),
  ('factor -> factor_int','factor',1,'p_factor','generic.py',226),
  ('factor_float -> signed_float','factor_float',1,'p_factor_float','generic.py',232),
  ('factor_float -> signed_float UINT signed_int','factor_float',3,'p_factor_float','generic.py',233),
  ('factor_float -> signed_float UINT power numeric_power','factor_float',4,'p_factor_float','generic.py',234),
  ('factor_int -> UINT','factor_int',1,'p_factor_int','generic.py',247),
  ('factor_int -> UINT signed_int','factor_int',2,'p_factor_int','generic.py',248),
  ('factor_int -> UINT power numeric_power','factor_int',3,'p_factor_int','generic.py',249),
  ('factor_int -> UINT UINT signed_int','factor_int',3,'p_factor_int','generic.py',250),
  ('factor_int -> UINT UINT power numeric_power','factor_int',4,'p_factor_int','generic.py',251),
  ('factor_fits -> UINT power OPEN_PAREN signed_int CLOSE_PAREN','factor_fits',5,'p_factor_fits','generic.py',269),
  ('factor_fits -> UINT power signed_int','factor_fits',3,'p_factor_fits','generic.py',270),
  ('factor_fits -> UINT SIGN UINT','factor_fits',3,'p_factor_fits','generic.py',271),
  ('factor_fits -> UINT OPEN_PAREN signed_int CLOSE_PAREN','factor_fits',4,'p_factor_fits','generic.py',272),
  ('product_of_units -> unit_expression product product_of_units','product_of_units',3,'p_product_of_units','generic.py',291),
  ('product_of_units -> unit_expression product_of_units','product_of_units',2,'p_product_of_units','generic.py',292),
  ('product_of_units -> unit_expression','product_of_units',1,'p_product_of_units','generic.py',293),
  ('unit_expression -> function','unit_expression',1,'p_unit_expression','generic.py',304),
  ('unit_expression -> unit_with_power','unit_expression',1,'p_unit_expression','generic.py',305),
  ('unit_expression -> OPEN_PAREN product_of_units CLOSE_PAREN','unit_expression',3,'p_unit_expression','generic.py',306),
  ('unit_with_power -> UNIT power numeric_power','unit_with_power',3,'p_unit_with_power','generic.py',315),
  ('unit_with_power -> UNIT numeric_power','unit_with_power',2,'p_unit_with_power','generic.py',316),
  ('unit_with_power -> UNIT','unit_with_power',1,'p_unit_with_power','generic.py',317),
  ('numeric_power -> sign UINT','numeric_power',2,'p_numeric_power','generic.py',328),
  ('numeric_power -> OPEN_PAREN paren_expr CLOSE_PAREN','numeric_power',3,'p_numeric_power','generic.py',329),
  ('paren_expr -> sign UINT','paren_expr',2,'p_paren_expr','generic.py',338),
  ('paren_expr -> signed_float','paren_expr',1,'p_paren_expr','generic.py',339),
  ('paren_expr -> frac','paren_expr',1,'p_paren_expr','generic.py',340),
  ('frac -> sign UINT division sign UINT','frac',5,'p_frac','generic.py',349),
  ('sign -> SIGN','sign',1,'p_sign','generic.py',355),
  ('sign -> <empty>','sign',0,'p_sign','generic.py',356),
  ('product -> STAR','product',1,'p_product','generic.py',365),
  ('product -> PERIOD','product',1,'p_product','generic.py',366),
  ('division -> SOLIDUS','division',1,'p_division','generic.py',372),
  ('power -> DOUBLE_STAR','power',1,'p_power','generic.py',378),
  ('power -> CARET','power',1,'p_power','generic.py',379),
  ('signed_int -> SIGN UINT','signed_int',2,'p_signed_int','generic.py',385),
  ('signed_float -> sign UINT','signed_float',2,'p_signed_float','generic.py',391),
  ('signed_float -> sign UFLOAT','signed_float',2,'p_signed_float','generic.py',392),
  ('function_name -> FUNCNAME','function_name',1,'p_function_name','generic.py',398),
  ('function -> function_name OPEN_PAREN main CLOSE_PAREN','function',4,'p_function','generic.py',404),
]
