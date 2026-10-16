# Publishing this repository on GitHub

The files in this ZIP are ready to become a new GitHub repository. They are **not yet published** and there is no live GitHub URL until you push them.

1. Sign in to https://github.com/new and create a repository named `local-ph-electrified-interfaces`. Choose public visibility only if you want the research code to be public. Do **not** initialize it with another README or `.gitignore`.
2. Extract the ZIP in your Downloads directory.
3. Run the commands below in a terminal from the extracted `wur-local-ph-github-english` folder (or whichever folder contains `README.md`).

```bash
cd ~/Downloads/wur-local-ph-github-english
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v

git init
git add .
git status --short
git commit -m "Publish local pH modeling and PHREEQC research portfolio"
git branch -M main
git remote add origin https://github.com/YOUR_USERNAME/local-ph-electrified-interfaces.git
git push -u origin main
```

Replace `YOUR_USERNAME` with your own GitHub handle. GitHub authentication may ask you to sign in via a credential manager or personal access token (never paste a token into a public chat or repository).

Before pushing, inspect `git status` and `git diff --cached --stat`. Do not publish CVs, motivation letters, private contact details or credentials. The repository includes `CITATION.cff`; review author attribution, and make a deliberate choice about licensing before inviting public reuse (`LICENSE_DECISION.md`).

Once published, the repository URL can be included in your CV. Do not describe it as a Wageningen research project or imply that university staff supervised it.
