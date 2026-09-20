# Notes on git from the command line for new users

These are notes for people unfamiliar with git and GitHub who need to clone the
repository — for beta testing, or to build from source.

GitHub no longer accepts username/password authentication over HTTPS, so you
need an SSH key. These notes go through it step by step.

> **This is written for a Linux environment.** On Windows use WSL2 (see
> [Building on Windows](building_on_windows.md)); on macOS the same commands
> work in Terminal or in a Multipass VM (see
> [Building on macOS](building_on_macos.md)).

---

## Clone using an SSH key

### 1. Generate an SSH key pair

Skip this if you already have one — check with `ls ~/.ssh/*.pub`.

```bash
ssh-keygen -t ed25519 -C "your_email@example.com"
# Press Enter to accept the default file location
# Add a passphrase for security
```

Replace `your_email@example.com` with your actual email address. The key pair
is saved to the default location: `~/.ssh/id_ed25519` (the **private** key,
never share it) and `~/.ssh/id_ed25519.pub` (the **public** key).

> Older guides use `-t rsa`, which produces `~/.ssh/id_rsa`. Ed25519 is the
> current recommendation — shorter, faster and at least as secure. If you have
> an existing `id_rsa` key it still works; just substitute the filename
> throughout.

### 2. Add the key to the SSH agent

The agent holds the decrypted key so you are not asked for the passphrase on
every operation. Start it first if it is not already running:

```bash
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_ed25519
```

The agent does not survive a reboot. If you get asked for your passphrase
repeatedly in later sessions, run these two lines again.

### 3. Copy your public key to GitHub

Print the **public** key:

```bash
cat ~/.ssh/id_ed25519.pub
```

Copy the entire line, including the `ssh-ed25519` prefix and the trailing
comment. Then:

1. Go to your GitHub profile settings.
2. Click **SSH and GPG keys**.
3. Click **New SSH key**.
4. Paste the key into the **Key** field.
5. Give it a title (e.g. "Work laptop" — something that tells you which machine
   it is, so you can revoke the right one later).
6. Click **Add SSH key**.

### 4. Check that it works

```bash
ssh -T git@github.com
```

Expect: `Hi <username>! You've successfully authenticated, but GitHub does not
provide shell access.` That message is a success, despite how it reads.

The first connection asks you to confirm GitHub's host fingerprint — answer
`yes`.

### 5. Clone the repository

```bash
git clone git@github.com:lsdtopotools/lsdtt3.git
```

More generally, replace the owner and repository name as needed:

```bash
git clone git@github.com:your_username/repository_name.git
```

Note the `git@github.com:` form — an SSH URL. The `https://github.com/...`
form is what triggers the password prompt that no longer works.

---

## Fixing a repository you already cloned over HTTPS

If you cloned before setting up SSH, the remote is still an HTTPS URL and will
keep asking for a password. Point it at the SSH URL instead:

```bash
git remote -v                                              # check what it is now
git remote set-url origin git@github.com:lsdtopotools/lsdtt3.git
git remote -v                                              # confirm the change
```

The subcommand is `git remote set-url`, not `git set-url`.

---

## Working on a branch

Never commit directly to `main`. Every change goes on its own branch, is
reviewed as a pull request, and is merged from there.

**The preferred workflow is: make a branch, work on it until you are satisfied,
open a pull request, and once it is merged delete the branch — then start the
next piece of work on a brand new branch.** Short-lived branches diverge less
from `main`, so the chance of a painful merge conflict stays low. A branch you
keep alive for weeks and keep adding unrelated work to is the main source of
conflicts in this repository.

### Make a branch

Start from an up-to-date `main`:

```bash
git switch main
git pull
git switch -c my-new-feature
```

`git switch -c` creates the branch and moves onto it. Name it after what the
work does — `fix-chi-nodata`, `add-hand-classifier` — not after yourself.

Then work as normal:

```bash
git add path/to/changed_file.cpp        # or: git add -A   to stage everything
git commit -m "Short description of what changed"
git push -u origin my-new-feature       # -u only needed on the first push
```

Commit as often as you like; commit messages on a branch are cheap.

### Open a pull request

Push the branch first, then either:

- **On GitHub.** After a push, GitHub shows a "Compare & pull request" button
  on the repository page. Click it, check the target is `main`, write a short
  description of what the change does and why, and create the pull request.
- **From the command line**, if you have the [GitHub CLI](https://cli.github.com/)
  installed:

  ```bash
  gh pr create --base main --fill
  gh pr view --web          # open it in a browser
  ```

Say in the description what you changed, what you tested, and anything a
reviewer should look at closely. If the pull request is not finished, open it
as a draft (`gh pr create --draft`) so nobody reviews it too early.

Push more commits to the same branch to update the pull request — it tracks the
branch, so there is no need to close and reopen it.

### Keep a branch in sync with `main`

If `main` moves on while you are working, bring those changes into your branch
so you find conflicts now rather than at merge time:

```bash
git switch main
git pull
git switch my-new-feature
git merge main
```

If git reports conflicts, work through them as described in
[Resolving merge conflicts](#resolving-merge-conflicts) below.

Rebasing (`git rebase main`) gives a tidier history but rewrites commits, which
means a force-push and confusion for anyone else looking at the branch. Use
`git merge` unless you know you want a rebase.

Sync early and sync often. The further behind `main` a branch drifts, the more
there is to reconcile.

### Resolving merge conflicts

A conflict happens when the same lines of a file changed on both `main` and
your branch, and git cannot decide which version to keep. It is not an error —
git is asking you a question. Nothing is lost, and the merge can always be
abandoned.

**1. See what is conflicted.**

```bash
git status                              # conflicts are listed under "Unmerged paths"
git diff --name-only --diff-filter=U    # just the file names
```

**2. Open each file and look for the markers.** git writes both versions into
the file:

```
<<<<<<< HEAD
float slope_threshold = 0.3;      // your branch's version
=======
float slope_threshold = 0.25;     // the version from main
>>>>>>> main
```

Everything between `<<<<<<< HEAD` and `=======` is what your branch has;
everything between `=======` and `>>>>>>>` is what you are merging in.

**3. Edit the file into the state you actually want**, and delete all three
marker lines. Sometimes you keep one side, sometimes the other, and often you
need to combine them — if someone else renamed a function on `main` and you
added a call to it, the correct answer is your new call using the new name.
Think about what the code should do; do not just pick a side to make the
markers go away.

**4. Mark it resolved and finish the merge.**

```bash
git add path/to/conflicted_file.cpp     # repeat for each file
git status                              # confirm nothing is left unmerged
git commit                              # completes the merge (a message is pre-filled)
git push
```

**5. Rebuild and run the tests before you push.** A file with the markers
removed still compiles as easily as one that is wrong. Conflicts are resolved
by hand, so the result is exactly as trustworthy as your rebuild:

```bash
cmake --build build -j
ctest --test-dir build -LE '^heavy$' --output-on-failure
```

Search the tree for stray markers too — a `<<<<<<<` left in a file is a
classic broken commit:

```bash
grep -rn '^<<<<<<< \|^>>>>>>> ' --include='*.cpp' --include='*.hpp' src/
```

#### If it goes wrong

Back out and start again — merging is not a one-way door:

```bash
git merge --abort                       # cancel the whole merge, branch untouched
git checkout --merge path/to/file.cpp   # restore the markers in one file you mangled
```

To take one side of a file wholesale, when you are certain that is right:

```bash
git checkout --ours   path/to/file.cpp  # keep your branch's version
git checkout --theirs path/to/file.cpp  # keep main's version
git add path/to/file.cpp
```

During a *merge*, "ours" is your branch and "theirs" is `main`. During a
*rebase* these are swapped, which is one more reason to prefer `git merge`.

A conflict you cannot judge is a conversation, not a puzzle. If the other side
of the conflict is somebody else's work and you cannot tell what they intended,
ask them before overwriting it.

#### Avoiding them in the first place

Most conflicts are a scheduling problem rather than a git problem. Merge `main`
into your branch often, keep branches short-lived, delete them after merging
(see above), and keep a branch to one piece of work so two people are less
likely to be editing the same lines at the same time.

### Delete the branch after merging

Once the pull request is merged, the branch has done its job. Delete it on
GitHub — there is a "Delete branch" button on the merged pull request — then
clean up locally:

```bash
git switch main
git pull                                # bring the merged work into local main
git branch -d my-new-feature            # delete the local branch
git fetch --prune                       # forget the deleted remote branch
```

`git branch -d` refuses to delete a branch whose work has not been merged,
which is a useful safety net. `git branch -D` forces it — only use that when
you are deliberately throwing the work away.

Start the next piece of work with a fresh `git switch -c` from the updated
`main`. Reusing an old branch for new work is what puts you back into conflict
territory.

---

## Executable bits on a Windows checkout

On a clone under `/mnt/c/` (NTFS via WSL2), `core.filemode` is effectively
`false`: NTFS has no executable bit, so `chmod +x` succeeds locally and git
sees nothing. `git status` stays clean and the mode never reaches the index.

This fails silently and for everyone else. A `debian/rules` or a shell script
committed as `100644` breaks the build for every person who clones, while
working perfectly on the machine it was committed from.

To set the bit in the repository itself:

```bash
git update-index --chmod=+x path/to/script.sh
git commit -m "Set executable bit on script.sh"
```

To check what the repository actually records, rather than what your working
tree shows:

```bash
git ls-tree HEAD debian/
```

`100755` is executable, `100644` is not. Do this for any new script, and
before releasing anything that runs `debian/rules` or `debian/smoke-test.sh`.

---

## See also

- [`building_from_source.md`](building_from_source.md) — what to do once you have cloned
- [`developer_onboarding.md`](developer_onboarding.md) — getting started as a contributor
