import pyperclip

def run():
    linetext = pyperclip.paste()
    if linetext == '\n':
        hashes = 80
    else:
        hashes = 80 - len(linetext)
    pyperclip.copy('#' * hashes)

run()
