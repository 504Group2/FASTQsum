html_template = '''<!doctype html>
<html>
    <head>
     <style>
        div{ Fastqsu,
        }       
     </style> 
</head>

<div class='p1' style="border: 1px solid blue">

#content1 This is my content
</div>

<div id='p1' style="border: 1px solid blue">
#content2
</div>

<div id='p1' style="border: 1px solid blue">
{scSum(csv)}
</div>

</html>
'''

with open('test3.html','w') as outf:
    outf.write(html_template.format(results=results))