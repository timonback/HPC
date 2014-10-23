#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014 Jakob Luettgau.
# Some rights reserved.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


import json
import os.path

###############################################################################
# settings
###############################################################################
raw_survey = {
    "id": "",
    "welcome": "HR 14/15 Rückmeldungen",
    "goodbye": "Danke!",
    "settings": {},
    "questions": [
        # id is optional
        {"id": "hrue-time", "question": "Bearbeitungszeit (in Stunden)?", "type": "numeric"},
        {"id": "hrue-diff", "question": "Schwierigkeit?\n(1 = zu leicht, 9 = zu schwer)", "type": "rate"},
        {"id": "hrue-salu", "question": "Lehrreich?\n(1 = wenig, 9 = sehr)", "type": "rate"},
        {"id": "hrue-arti", "question": "Verständlichkeit?\n(1 = großteils unklar, 9 = verständlich)", "type": "rate"},
        {"id": "hrue-comm", "question": """Kommentar 
   z.B. 
   - Zeitaufwand zu groß?
   - Wo seid ihr nicht weiter gekommen?
   - Was war zu schwer oder unverständlich?
   - Was könnten wir verbessern?
   - Was hat euch besonders gut gefallen?"""},
    ],
}

output = "./feedback.txt"
###############################################################################

# helpers
class Survey(object):
    def __init__(self, survey={}, outfile="./feedback.txt"):
        self.survey = survey
        self.outfile = outfile

    def rate(self, question):
        if not "hint" in question:
            question["hint"] = "(1-9) "
        return self.prompt(question, options={"range": range(1,10)})

    def numeric(self, question):
        if not "hint" in question:
            question["hint"] = "(numeric) "
        return self.prompt(question, options={"allow": "int"})

    def choice(self, question, choices=[]):
        if not "hint" in question:
            question["hint"] = choices
        return self.prompt(question)

    def yesno(self, question):
        return self.choice(question, ["yes", "no"])

    def text(self, question):
        return self.prompt(question)
        pass

    # interactive
    def begin(self):
        # say hello
        if "welcome" in self.survey:
            print self.survey["welcome"], "\n"
        # ask questions
        for q in self.survey["questions"]:
            # print question
            print "Q:", q["question"]
            # ensure a question type is set
            if not "type" in q:
                q["type"] = "text"
            # handle type and prompt user
            if q["type"] == "rate":
                q["answer"] = self.rate(q)
            elif q["type"] == "choice":
                q["answer"] = self.choice(q)
            elif q["type"] == "yesno":
                q["answer"] = self.yesno(q)
            elif q["type"] == "numeric":
                q["answer"] = self.numeric(q)
            else:
                # fallback to text
                q["answer"] = self.text(q)
            # leave some room
            print ""

        # print goodbye
        if "goodbye" in self.survey:
            print self.survey["goodbye"]


    def prompt(self, question, options={}):
        """ Handle different kind of prompts
            @param  mode
            @return sanitized user input
        """
        
        if "hint" in question:
            hint = question["hint"]
        else:
            hint = ""

        if "print_question" in options:
            if "question" in question:
                print "Q:", question["question"]

        print "A:", 
        input = raw_input(hint)

        # sanitize choice
        if isinstance(hint, list):
            if input in hint:
                # all good
                return input
            else:
                return self.prompt(question)

        # sanitze range
        if "range" in options:
            try:
                input = int(input)
            except:
                pass

            if input in options["range"]:
                # all good
                return input
            else:
                return self.prompt(question, options)

        # sanitze numeric input
        if "allow" in options:
            if options["allow"] == "int":
                try:
                    input = int(input)
                    return input
                except:
                    return self.prompt(question, options)

        # discourage short answers
        if not "optional" in options:
            if len(input) < 3:
                print "Seems to be a little bit short. Try again?"
                if self.yesno({}) == "yes":
                    options["print_question"] = True
                    return self.prompt(question, options)

        return input

    # dump
    def save(self):
        buf = json.dumps(self.survey, ensure_ascii=False)
        outfile = self.outfile
        # overwrite protection
        counter = 0
        while os.path.isfile(outfile):
            counter += 1
            outfile = self.outfile + "." + str(counter)
        f = open(outfile, 'w')
        f.write(buf)
        f.close()
        pass


# logic
def main():
    s = Survey(raw_survey, output)
    s.begin()
    s.save()
    pass

if __name__ == '__main__':
    main()
