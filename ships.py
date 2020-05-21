from retrograde import *

# shipyard
dinghy = [
	Module('generator',[0,0]),
	Module('engine 10',[0,15]),
	Module('RCS',[0,-10]) ]

lothar = [
	Module('generator',[0,0]),
	Module('engine 10',[-13,8],0.6/math.pi),
	Module('engine 10',[13,8],-0.6/math.pi),
	Module('RCS',[-13,-10]),
	Module('RCS',[13,-10]) ,
	Module('cannon 10',[0,-10]) ]

boldang = [
	Module('spar 10',[0,-100],(0.5* math.pi)),
	Module('box 10',[0,0])]

bigmolly = [
	Module('box 100',[0,0]),
	Module('spar 100',[1000,0],0.5 * math.pi),
	Module('box 100',[-1000,0]),
	Module('box 100',[2000,0]),
	Module('box 100',[-2000,0]),
	Module('box 100',[3000,0])]

derelict_hyperunit = [
	Module('hyperdrive 10',[0,0])]

ida_frigate = [
	Module('generator',[0,50]),
	Module('engine 10',[-13,58],0.6/math.pi),
	Module('engine 10',[13,58],-0.6/math.pi),
	Module('RCS',[-13,40]),
	Module('RCS',[13,40]) ,
	Module('box 10',[0,10]),
	Module('box 10',[0,-40]),
	Module('hyperdrive 10',[0,-75]),
	Module('RCS',[-13,-70]),
	Module('RCS',[13,-70]) ]
