s1 = set(
[str(sorted([u'102', u'101', u'104'])),
str(sorted([u'102', u'101', u'118'])),
str(sorted([u'102', u'118', u'104'])),
str(sorted([u'107', u'109', u'106'])),
str(sorted([u'107', u'111', u'106'])),
str(sorted([u'107', u'114', u'106'])),
str(sorted([u'107', u'114', u'109'])),
str(sorted([u'107', u'114', u'111'])),
str(sorted([u'107', u'114', u'118'])),
str(sorted([u'109', u'102', u'118'])),
str(sorted([u'109', u'106', u'118'])),
str(sorted([u'109', u'111', u'106'])),
str(sorted([u'114', u'102', u'118'])),
str(sorted([u'114', u'106', u'118'])),
str(sorted([u'114', u'109', u'106'])),
str(sorted([u'114', u'109', u'118'])),
str(sorted([u'114', u'111', u'106'])),
str(sorted([u'114', u'122', u'118'])),
str(sorted([u'120', u'102', u'101'])),
str(sorted([u'120', u'102', u'104'])),
str(sorted([u'120', u'102', u'118'])),
str(sorted([u'120', u'122', u'102'])),
str(sorted([u'122', u'102', u'101'])),
str(sorted([u'122', u'102', u'104'])),
str(sorted([u'122', u'102', u'118'])),
str(sorted([u'122', u'109', u'118']))
]
)

s2 = set([str(sorted([u'101', u'102', u'118'])),
str(sorted([u'101', u'102', u'120'])),
str(sorted([u'101', u'102', u'122'])),
str(sorted([u'102', u'104', u'118'])),
str(sorted([u'102', u'104', u'120'])),
str(sorted([u'102', u'104', u'122'])),
str(sorted([u'102', u'118', u'120'])),
str(sorted([u'102', u'118', u'122'])),
str(sorted([u'102', u'120', u'122'])),
str(sorted([u'106', u'107', u'109'])),
str(sorted([u'106', u'107', u'111'])),
str(sorted([u'106', u'107', u'114'])),
str(sorted([u'106', u'109', u'111'])),
str(sorted([u'106', u'109', u'114'])),
str(sorted([u'106', u'109', u'118'])),
str(sorted([u'106', u'111', u'114'])),
str(sorted([u'106', u'114', u'118'])),
str(sorted([u'107', u'111', u'114'])),
str(sorted([u'107', u'114', u'118'])),
str(sorted([u'109', u'114', u'118']))])

print s1 - s2
print len(s1 - s2)



